/**
 * License: GPL
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 2
 * as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
package org.janelia.saalfeldlab;

import java.util.Collection;

import org.janelia.saalfeldlab.n5.GzipCompression;
import org.janelia.saalfeldlab.n5.N5Writer;
import org.janelia.saalfeldlab.n5.hdf5.N5HDF5Writer;
import org.janelia.saalfeldlab.n5.imglib2.N5Utils;

import bdv.bigcat.annotation.Annotation;
import bdv.bigcat.annotation.Annotations;
import bdv.bigcat.annotation.AnnotationsHdf5Store;
import bdv.bigcat.annotation.PostSynapticSite;
import bdv.bigcat.annotation.PreSynapticSite;
import bdv.util.LocalIdService;
import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.Interval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealLocalizable;
import net.imglib2.RealPoint;
import net.imglib2.img.cell.CellImgFactory;
import net.imglib2.interpolation.randomaccess.NLinearInterpolatorFactory;
import net.imglib2.realtransform.AffineGet;
import net.imglib2.realtransform.AffineRealRandomAccessible;
import net.imglib2.realtransform.RealViews;
import net.imglib2.realtransform.Scale3D;
import net.imglib2.type.numeric.integer.UnsignedByteType;
import net.imglib2.type.numeric.integer.UnsignedLongType;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;
import picocli.CommandLine;
import picocli.CommandLine.Option;

/**
 *
 *
 * @author Stephan Saalfeld &lt;saalfelds@janelia.hhmi.org&gt;
 */
public class Centrosomes {

	public static class Parameters {

		@Option(names = { "--infile", "-i" }, description = "input CREMI-format HDF5 file name (default == --infile")
		public String inFile = null;
	}

	private static final RealPoint getPosition(
			final Annotation a,
			final Interval interval,
			final double[] offset,
			final double[] resolution) {

		final int n = interval.numDimensions();
		final RealPoint aPosition = a.getPosition();
		final RealPoint position = new RealPoint(n);
		for (int d = 0; d < n; ++d)
			position.setPosition(Math.max(0, Math.min(interval.max(d), Math.round((aPosition.getDoublePosition(d) - offset[d]) / resolution[d]))), d);
		return position;
	}

	static final private Interval boundingBox(
			final RealLocalizable a,
			final RealLocalizable b,
			final double radius,
			final double[] offset,
			final double[] resolution,
			final Interval src) {

		final int n = a.numDimensions();
		final long[] min = new long[n];
		final long[] max = new long[n];
		for (int d = 0; d < n; ++d) {
			final double aPos = a.getDoublePosition(d) - offset[d];
			final double bPos = b.getDoublePosition(d) - offset[d];
			min[d] = Math.round(Math.max(src.min(d), (Math.min(aPos, bPos) - radius) / resolution[d]));
			max[d] = Math.round(Math.min(src.max(d), (Math.max(aPos, bPos) + radius) / resolution[d]));
		}
		return new FinalInterval(min, max);
	}

	private static final double squareLength(final RealLocalizable a) {

		double l = 0;
		for (int d = 0; d < a.numDimensions(); ++d)
			l += a.getDoublePosition(d) * a.getDoublePosition(d);
		return l;
	}

	private static final void norm(final RealPoint a) {

		final double l = Math.sqrt(squareLength(a));

		for (int d = 0; d < a.numDimensions(); ++d)
			a.setPosition(a.getDoublePosition(d) / l, d);
	}

	private static final double dot(
			final RealLocalizable a,
			final RealLocalizable b) {

		double l = 0;
		for (int d = 0; d < a.numDimensions(); ++d)
			l += a.getDoublePosition(d) * b.getDoublePosition(d);

		return l;
	}

	private static final void cross(
			final RealLocalizable a,
			final RealLocalizable b,
			final RealPoint c) {

		final double x = a.getDoublePosition(1) * b.getDoublePosition(2) - a.getDoublePosition(2) * b.getDoublePosition(1);
		final double y = a.getDoublePosition(2) * b.getDoublePosition(0) - a.getDoublePosition(0) * b.getDoublePosition(2);
		final double z = a.getDoublePosition(0) * b.getDoublePosition(1) - a.getDoublePosition(1) * b.getDoublePosition(0);

		c.setPosition(x, 0);
		c.setPosition(y, 1);
		c.setPosition(z, 2);
	}

	private static final void subtract(
			final RealLocalizable a,
			final RealLocalizable b,
			final RealPoint c) {

		for (int d = 0; d < a.numDimensions(); ++d)
			c.setPosition(a.getDoublePosition(d) - b.getDoublePosition(d), d);
	}

	private static double squareDistanceToLine(
			final RealLocalizable p,
			final RealLocalizable a,
			final RealLocalizable b,
			final RealPoint c) {

		subtract(p, a, c);
		cross(b, c, c);
		return squareLength(c) / squareLength(b);
	}

	private static boolean inCylinder(
			final RealLocalizable p,
			final RealLocalizable a,
			final RealLocalizable b,
			final RealPoint c,
			final RealPoint d,
			final double squareRadius) {

		subtract(b, a, d);
		subtract(p, a, c);
		final double squareLengthB = squareLength(d);
		final double t = dot(c, d) / squareLengthB;

		if (t >= 0 && t <= 1) {
			cross(d, c, c);
			return squareLength(c) / squareLengthB <= squareRadius;
		}
		return false;
	}

	private static void paintMicrotubule(
			final RandomAccessibleInterval<UnsignedLongType> src,
			final RealPoint a,
			final RealPoint b,
			final double outerRadius,
			final double innerRadius,
			final double[] offset,
			final double[] resolution,
			final long fg,
			final long bg) {

		final IntervalView<UnsignedLongType> box = Views.interval(src, boundingBox(a, b, outerRadius, offset, resolution, src));
		final RealPoint p = new RealPoint(src.numDimensions());
		final RealPoint c = new RealPoint(src.numDimensions());
		final RealPoint d = new RealPoint(src.numDimensions());
		final double sor = outerRadius * outerRadius;
		final double sir = innerRadius * innerRadius;
		final Cursor<UnsignedLongType> cursor = box.cursor();
		while (cursor.hasNext()) {
			cursor.fwd();
			for (int i = 0; i < src.numDimensions(); ++i)
				p.setPosition(cursor.getDoublePosition(i) * resolution[i] + offset[i], i);
			if (inCylinder(p, a, b, c, d, sor))
			if (true)
				cursor.get().set(fg);
		}
		cursor.reset();
		while (cursor.hasNext()) {
			cursor.fwd();
			for (int i = 0; i < src.numDimensions(); ++i)
				p.setPosition(cursor.getDoublePosition(i) * resolution[i] + offset[i], i);
			if (inCylinder(p, a, b, c, d, sir))
				cursor.get().set(bg);
		}
	}

	private static void reverse(final double[] a) {

		for (int i = a.length / 2, j = a.length - 1 - i; i >= 0; --i, ++j) {
			final double ai = a[i];
			a[i] = a[j];
			a[j] = ai;
		}
	}

	private static void paintMicrotubule(
			final Annotation a,
			final RandomAccessibleInterval<UnsignedLongType> centrosomesVolume,
			final double[] fLabelsOffset,
			final double[] labelsResolution) {

		final PreSynapticSite b = ((PostSynapticSite)a).getPartner();
		final RealPoint posA = a.getPosition();
		final RealPoint posB = b.getPosition();

		paintMicrotubule(centrosomesVolume, posA, posB, 21, 18, fLabelsOffset, labelsResolution, 255, 0);

		System.out.println(String.format(
				"%s -> %s: (%f, %f, %f) -> (%f, %f, %f)",
				b.getComment(),
				a.getComment(),
				posB.getDoublePosition(0),
				posB.getDoublePosition(1),
				posB.getDoublePosition(2),
				posA.getDoublePosition(0),
				posA.getDoublePosition(1),
				posA.getDoublePosition(2)));
	}




	/**
	 * @param args
	 * @throws Exception
	 */
	public static void main(final String... args) throws Exception {

//		new ImageJ();

		final Parameters params = new Parameters();
		try {
			CommandLine.populateCommand(params, args);
		} catch (final RuntimeException e) {
			CommandLine.usage(params, System.err);
			return;
		}

//		System.out.println("Opening input " + params.inFile);

		final IHDF5Writer hdf5Writer = HDF5Factory.open(params.inFile);
		final N5Writer n5 = new N5HDF5Writer(hdf5Writer, new int[]{64, 64, 64});

		//final N5HDF5Reader n5 = new N5HDF5Reader(params.inFile, new int[]{64, 64, 64});

		final RandomAccessibleInterval<UnsignedLongType> labelsSource = N5Utils.open(n5, "/volumes/labels/gt");
		final double[] labelsResolution = n5.getAttribute("/volumes/labels/gt", "resolution", double[].class);
		reverse(labelsResolution);
		final double[] labelsOffset = n5.getAttribute("/volumes/labels/gt", "offset", double[].class);
		reverse(labelsOffset);
		final double[] fLabelsOffset =  (labelsOffset == null) ? new double[3] : labelsOffset;

		System.out.println(String.format(
				"label offset : (%f, %f, %f)",
				fLabelsOffset[0],
				fLabelsOffset[1],
				fLabelsOffset[2]));

		System.out.println(String.format(
				"label resolution : (%f, %f, %f)",
				labelsResolution[0],
				labelsResolution[1],
				labelsResolution[2]));

		final RandomAccessibleInterval<UnsignedLongType> centrosomesVolume =
				new CellImgFactory<>(new UnsignedLongType(), new int[] {64, 64, 64}).create(labelsSource);

		/* annotations */
		final AnnotationsHdf5Store annotationsStore = new AnnotationsHdf5Store(params.inFile, new LocalIdService());
		final Annotations annotations = annotationsStore.read();

		final Collection<Annotation> annotationsCollection = annotations.getAnnotations();

		annotationsCollection.forEach(
				a -> {
					if (a instanceof PostSynapticSite && a.getComment().equalsIgnoreCase("c"))
						paintMicrotubule(a, centrosomesVolume, fLabelsOffset, labelsResolution);
				});
		annotationsCollection.forEach(
				a -> {
					if (a instanceof PostSynapticSite && a.getComment().equalsIgnoreCase("b"))
						paintMicrotubule(a, centrosomesVolume, fLabelsOffset, labelsResolution);
				});
		annotationsCollection.forEach(
				a -> {
					if (a instanceof PostSynapticSite && a.getComment().equalsIgnoreCase("a"))
						paintMicrotubule(a, centrosomesVolume, fLabelsOffset, labelsResolution);
				});

		RandomAccessibleInterval<UnsignedByteType> raw = N5Utils.open(n5, "/volumes/raw");
		final double[] rawResolution = n5.getAttribute("/volumes/raw", "resolution", double[].class);
		reverse(rawResolution);
		raw = Views.translate(
				raw,
				-Math.round(fLabelsOffset[0] / rawResolution[0]),
				-Math.round(fLabelsOffset[1] / rawResolution[1]),
				-Math.round(fLabelsOffset[2] / rawResolution[2]));
		final Scale3D rawScale = new Scale3D(
				rawResolution[0] / labelsResolution[0],
				rawResolution[1] / labelsResolution[1],
				rawResolution[2] / labelsResolution[2]);
		final AffineRealRandomAccessible<UnsignedByteType, AffineGet> rawScaled = RealViews.affineReal(Views.interpolate(Views.extendZero(raw), new NLinearInterpolatorFactory<>()), rawScale);


//		final BdvSource bdv = BdvFunctions.show(rawScaled, labelsSource, "raw");
//		final BdvOptions options = BdvOptions.options().addTo(bdv);
//		BdvFunctions.show(labelsSource, "labels", options);
//		BdvFunctions.show(centrosomesVolume, "centrosomes", options);

		N5Utils.save(centrosomesVolume, n5, "/volumes/labels/centrosomes", new int[] {64, 64, 64}, new GzipCompression());
		reverse(fLabelsOffset);
		reverse(labelsResolution);
		n5.setAttribute("/volumes/labels/centrosomes", "offset", fLabelsOffset);
		n5.setAttribute("/volumes/labels/centrosomes", "resolution", labelsResolution);

		hdf5Writer.close();

//		System.out.println("Done.");
	}

}

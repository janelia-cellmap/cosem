/**
 * THE CRAPL v0 BETA 1
 */
package org.janelia.saalfeldlab;

import java.io.IOException;
import java.util.Arrays;
import java.util.Optional;
import java.util.concurrent.Callable;

import org.janelia.saalfeldlab.n5.GzipCompression;
import org.janelia.saalfeldlab.n5.N5Reader;
import org.janelia.saalfeldlab.n5.N5Writer;
import org.janelia.saalfeldlab.n5.ij.N5Factory;
import org.janelia.saalfeldlab.n5.imglib2.N5Utils;

import bdv.util.BdvOptions;
import bdv.util.volatiles.SharedQueue;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.interpolation.randomaccess.NLinearInterpolatorFactory;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.realtransform.RealViews;
import net.imglib2.realtransform.Scale3D;
import net.imglib2.realtransform.Translation3D;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.util.Intervals;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;
import picocli.CommandLine;
import picocli.CommandLine.Option;

public class CropBigCat implements Callable<Void> {

	@Option(
			names = {"--rawn5url"},
			required = true,
			description = "N5 URL to raw volume, e.g. '--rawn5url /groups/cosem/cosem/annotations/amira/hemibrain_crop153_roi1/test712x712x712+124+219+54.h5'")
	private String rawN5Url = null;

//	@Option(
//			names = {"--labelsn5url"},
//			required = false,
//			description = "N5 URL to labels volume, e.g. '--labelsn5url /groups/cosem/cosem/annotations/amira/hemibrain_crop153_roi1/test712x712x712+124+219+54.h5'")
//	private String labelsN5Url = null;

	@Option(
			names = {"--outn5url"},
			required = false,
			description = "N5 URL to export BigCAT project, e.g. '--outn5url /groups/cosem/cosem/data/hemi-brain_8x8x8nm/hemibrain_crop153_roi1_reoriented.h5'")
	private String outN5Url = null;

	@Option(
			names = {"--rawres"},
			required = false,
			split = ",",
			description = "Raw resolution as a comma separated list of doubles, if present overrides resolution attribute in raw dataset, e.g. --rawres 4,4,4")
	final double[] rawResolution = null;

	@Option(
			names = {"--rawoffset"},
			required = false,
			split = ",",
			description = "Raw offset as a comma separated list of doubles in world spacings, if present overrides offset attribute in raw dataset, e.g. --rawoffset 0,0,0")
	final double[] rawOffset = null;

	@Option(
			names = {"-r", "--rawdataset"},
			required = true,
			description = "raw dataset, e.g. '-r /volumes/raw'")
	private String rawPath = null;

//	@Option(
//			names = {"-l", "--labelsdataset"},
//			required = true,
//			description = "labels dataset, e.g. '-l /volumes/labels/gt'")
//	private String labelsPath = null;
//
//	@Option(
//			names = {"--labelsres"},
//			required = false,
//			split = ",",
//			description = "Labels resolution as a comma separated list of doubles, if present overrides resolution attribute in labels dataset, e.g. --labelsres 4,4,4")
//	final double[] labelsResolution = null;
//
//	@Option(
//			names = {"--labelsoffset"},
//			required = false,
//			split = ",",
//			description = "Labels offset as a comma separated list of doubles in world spacings, if present overrides offset attribute in labels dataset, e.g. --labelsoffset 0,0,0")
//	final double[] labelsOffset = null;

	@Option(
			names = {"-b", "--bookmark"},
			required = true,
			split = ",",
			description = "BigDataViewer bookmark as a comma separated list of doubles, e.g. -b 4,4,4,4,4,4,4,4,4,4,4,4")
	final double[] bookmark = null;

	@Option(
			names = {"--labelssize"},
			required = true,
			split = ",",
			description = "export size as a comma separated list of longs in pixels, e.g. --labelssize 1024,1024,1024")
	final long[] labelsSize = null;

	@Option(
			names = {"--rawpadding"},
			required = true,
			split = ",",
			description = "raw padding size as a comma separated list of longs in pixels, e.g. --rawpadding 256,256,256")
	final long[] rawPadding = null;



	/**
	 * Start the tool.  We ignore the exit code returned by
	 * {@link CommandLine#execute(String...)} but this can be useful in other
	 * applications.
	 *
	 * @param args
	 */
	public static void main(final String... args) {

		new CommandLine(new CropBigCat()).execute(args);
	}

	/**
	 * The real implementation.  We use {@link Callable Callable<Void>} instead
	 * of {@link Runnable} because {@link Runnable#run()} cannot throw
	 * {@link Exception Exceptions}.
	 *
	 * Since we would like to use some type parameters, we have to delegate to
	 * a method that was not declared in an interface without such parameters.
	 *
	 * @throws Exception
	 */
	@Override
	public Void call() throws Exception {

		run();
		return null;
	}

	private static double axisNorm(final AffineTransform3D affine, final int d) {

		final double x = affine.get(0, d);
		final double y = affine.get(1, d);
		final double z = affine.get(2, d);

		return Math.sqrt(x * x + y * y + z * z);
	}

	private final <R extends NativeType<R> & RealType<R>, L extends NativeType<L> & RealType<L>> void run() throws IOException {

		final SharedQueue queue = new SharedQueue(Math.max(1, Runtime.getRuntime().availableProcessors() - 1));

		final N5Reader n5Raw = new N5Factory().openReader(rawN5Url);
		final RandomAccessibleInterval<R> raw = N5Utils.openVolatile(n5Raw, rawPath);

//		final N5Reader n5Labels = new N5Factory().openReader(Optional.ofNullable(labelsN5Url).orElse(rawN5Url));
//		final RandomAccessibleInterval<L> labels = N5Utils.openVolatile(n5Labels, labelsPath);

		final double[] rawResolution =
				Optional.ofNullable(this.rawResolution)
				.orElse(
						Optional.ofNullable(n5Raw.getAttribute(rawPath, "resolution", double[].class))
						.orElse(new double[]{1, 1, 1}));
		final double[] rawOffset =
				Optional.ofNullable(this.rawOffset)
				.orElse(
						Optional.ofNullable(n5Raw.getAttribute(rawPath, "offset", double[].class))
						.orElse(new double[]{0, 0, 0}));
//		final double[] labelsResolution =
//				Optional.ofNullable(this.labelsResolution)
//				.orElse(
//						Optional.ofNullable(n5Labels.getAttribute(labelsPath, "resolution", double[].class))
//						.orElse(new double[]{1, 1, 1}));
//		final double[] labelsOffset =
//				Optional.ofNullable(this.labelsOffset)
//				.orElse(
//						Optional.ofNullable(n5Labels.getAttribute(labelsPath, "offset", double[].class))
//						.orElse(new double[]{0, 0, 0}));

		System.out.println("rawResolution = " + Arrays.toString(rawResolution));
		System.out.println("rawOffset = " + Arrays.toString(rawOffset));
//		System.out.println("labelsResolution = " + Arrays.toString(labelsResolution));
//		System.out.println("labelsOffset = " + Arrays.toString(labelsOffset));

		final Scale3D rawScale = new Scale3D(rawResolution);
		final Translation3D rawTranslation = new Translation3D(rawOffset);

		final Scale3D labelsScale = new Scale3D(new double[] {4, 4, 4});
//		final Scale3D labelsScale = new Scale3D(labelsResolution);
//		final Translation3D labelsTranslation = new Translation3D(labelsOffset);

		final AffineTransform3D bookmarkAffine = new AffineTransform3D();
		bookmarkAffine.set(bookmark);

		final double norm0 = axisNorm(bookmarkAffine, 0);
		final double norm1 = axisNorm(bookmarkAffine, 1);
		final double norm2 = axisNorm(bookmarkAffine, 2);

		System.out.println(norm0);
		System.out.println(norm1);
		System.out.println(norm2);

		final Scale3D scale = new Scale3D(norm0, norm1, norm2);

		final Translation3D labelsCenterOffset = new Translation3D(
				0.5 * labelsSize[0],
				0.5 * labelsSize[1],
				0.5 * labelsSize[2]);

		final Translation3D rawCenterOffset = new Translation3D(
				0.5 * labelsSize[0] + rawPadding[0],
				0.5 * labelsSize[1] + rawPadding[1],
				0.5 * labelsSize[2] + rawPadding[2]);

		bookmarkAffine.preConcatenate(scale.inverse());
		bookmarkAffine.preConcatenate(labelsScale.inverse());

		final AffineTransform3D rawAffine = new AffineTransform3D();
		rawAffine.preConcatenate(rawScale);
		rawAffine.preConcatenate(rawTranslation);
//		rawAffine.preConcatenate(rawScale.inverse());
		rawAffine.preConcatenate(bookmarkAffine);
		rawAffine.preConcatenate(rawCenterOffset);

//		final AffineTransform3D labelsAffine = new AffineTransform3D();
//		labelsAffine.preConcatenate(labelsScale);
//		labelsAffine.preConcatenate(labelsTranslation);
//		labelsAffine.preConcatenate(bookmarkAffine);
//		labelsAffine.preConcatenate(labelsCenterOffset);

		final BdvOptions options = BdvOptions.options();

		final IntervalView<R> rawDeformed = Views.interval(
				RealViews.affine(
					Views.interpolate(
							Views.extendZero(raw),
							new NLinearInterpolatorFactory<>()),
//							new NearestNeighborInterpolatorFactory<>()),
					rawAffine),
				Intervals.createMinSize(
						0, 0, 0,
						labelsSize[0] + rawPadding[0] + rawPadding[0],
						labelsSize[1] + rawPadding[1] + rawPadding[1],
						labelsSize[2] + rawPadding[2] + rawPadding[2]));

//		final IntervalView<L> labelsDeformed = Views.interval(
//				RealViews.affine(
//					Views.interpolate(
//							Views.extendZero(labels),
//	//						new NLinearInterpolatorFactory<>()),
//							new NearestNeighborInterpolatorFactory<>()),
//					labelsAffine),
//				Intervals.createMinSize(
//						0, 0, 0,
//						labelsSize[0],
//						labelsSize[1],
//						labelsSize[2]));

		final N5Writer n5Writer = new N5Factory().openWriter(outN5Url);
		N5Utils.save(
				rawDeformed,
				n5Writer,
				"/volumes/raw",
				new int[] {64, 64, 64},
				new GzipCompression());
//		N5Utils.save(
//				labelsDeformed,
//				n5Writer,
//				labelsPath,
//				new int[] {64, 64, 64},
//				new GzipCompression());
		n5Writer.setAttribute("/volumes/raw", "resolution", rawResolution);
		n5Writer.setAttribute("/volumes/raw", "offset", new double[3]);
//		n5Writer.setAttribute(labelsPath, "resolution", rawResolution);
//		n5Writer.setAttribute(
//				labelsPath,
//				"offset",
//				new double[] {
//						rawPadding[0] * rawResolution[0],
//						rawPadding[1] * rawResolution[1],
//						rawPadding[2] * rawResolution[2]});

		n5Writer.setAttribute("/", "rawN5Url", rawN5Url);
//		if (labelsN5Url != null) n5Writer.setAttribute("/", "labelsN5Url", labelsN5Url);
		if (this.rawResolution != null) n5Writer.setAttribute("/", "rawResolution", this.rawResolution);
		if (this.rawOffset != null) n5Writer.setAttribute("/", "rawOffset", this.rawOffset);
		n5Writer.setAttribute("/", "rawPath", rawPath);
//		n5Writer.setAttribute("/", "labelsPath", labelsPath);
//		if (this.labelsResolution != null) n5Writer.setAttribute("/", "labelsResolution", this.labelsResolution);
//		if (this.labelsOffset != null) n5Writer.setAttribute("/", "labelsOffset", this.labelsOffset);
		n5Writer.setAttribute("/", "bookmark", bookmark);
		n5Writer.setAttribute("/", "labelsSize", labelsSize);
		n5Writer.setAttribute("/", "rawPadding", rawPadding);

		n5Writer.close();

//		final BdvStackSource<R> bdv = BdvFunctions.show(
////				rawDeformed,
//				RealViews.affineReal(
//						Views.interpolate(
//								Views.extendZero(raw),
//		//						new NLinearInterpolatorFactory<>()),
//								new NearestNeighborInterpolatorFactory<>()),
//						rawAffine),
//				Intervals.createMinMax(-1, -1, -1, 1, 1, 1),
//				rawPath);
//
//		BdvFunctions.show(
////				labelsDeformed,
//				RealViews.affineReal(
//						Views.interpolate(
//								Views.extendZero(labels),
////								new NLinearInterpolatorFactory<>()),
//								new NearestNeighborInterpolatorFactory<>()),
//						labelsAffine),
//				Intervals.createMinMax(-1, -1, -1, 1, 1, 1),
//				labelsPath,
//				options.addTo(bdv));

//		final BdvStackSource<R> bdv1 = BdvFunctions.show(
//				Views.offsetInterval(
//						raw,
//						Intervals.createMinSize(
//								Math.round(-rawOffset[0] / rawResolution[0]),
//								Math.round(-rawOffset[1] / rawResolution[1]),
//								Math.round(-rawOffset[2] / rawResolution[2]),
//								1024,
//								1024,
//								1024)),
//				"");

//		final BdvStackSource<R> bdv1 = BdvFunctions.show(
//				N5Utils.<R>open(n5Labels, "/volumes/raw"),
//				"hdf5");
//
//		BdvFunctions.show(
////				rawDeformed,
//				RealViews.affineReal(
//						Views.interpolate(
//								Views.extendZero(raw),
//		//						new NLinearInterpolatorFactory<>()),
//								new NearestNeighborInterpolatorFactory<>()),
//						rawAffine),
//				Intervals.createMinMax(-1, -1, -1, 1, 1, 1),
//				rawPath,
//				options.addTo(bdv1));

		System.out.println("Done.");

	}
}

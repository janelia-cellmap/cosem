/**
 *
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
import bdv.util.LocalIdService;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealPoint;
import net.imglib2.img.cell.CellImgFactory;
import net.imglib2.type.numeric.integer.UnsignedLongType;
import picocli.CommandLine;
import picocli.CommandLine.Option;

/**
 * Export all presynaptic partner annotations as a single pixel.
 *
 * @author Stephan Saalfeld &lt;saalfelds@janelia.hhmi.org&gt;
 */
public class Ribosomes {

	public static class Parameters {

		@Option(names = { "--infile", "-i" }, description = "input CREMI-format HDF5 file name (default == --infile")
		public String inFile = null;
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

		final N5Writer n5 = new N5HDF5Writer(params.inFile, new int[]{64, 64, 64});
		final RandomAccessibleInterval<UnsignedLongType> labelsSource = N5Utils.open(n5, "/volumes/labels/gt");
		final double[] labelsResolution = n5.getAttribute("/volumes/labels/gt", "resolution", double[].class);
		final double[] labelsOffset = n5.getAttribute("/volumes/labels/gt", "offset", double[].class);
		final double[] fLabelsOffset =  (labelsOffset == null) ? new double[3] : labelsOffset;

		final RandomAccessibleInterval<UnsignedLongType> ribosomesVolume =
				new CellImgFactory<>(new UnsignedLongType(), new int[] {64, 64, 64}).create(labelsSource);
		final RandomAccess<UnsignedLongType> ribosomesVolumeAccess = ribosomesVolume.randomAccess();

		/* annotations */
		final AnnotationsHdf5Store annotationsStore = new AnnotationsHdf5Store(params.inFile, new LocalIdService());
		final Annotations annotations = annotationsStore.read();

		final Collection<Annotation> annotationsCollection = annotations.getAnnotations();

		annotationsCollection.forEach(
				a -> {
					final RealPoint position = a.getPosition();
					ribosomesVolumeAccess.setPosition(
							Math.max(0, Math.min(ribosomesVolume.max(0), Math.round((position.getDoublePosition(0) - fLabelsOffset[0]) / labelsResolution[0]))), 0);
					ribosomesVolumeAccess.setPosition(
							Math.max(0, Math.min(ribosomesVolume.max(1), Math.round((position.getDoublePosition(1) - fLabelsOffset[1]) / labelsResolution[1]))), 1);
					ribosomesVolumeAccess.setPosition(
							Math.max(0, Math.min(ribosomesVolume.max(2), Math.round((position.getDoublePosition(2) - fLabelsOffset[2]) / labelsResolution[2]))), 2);
					ribosomesVolumeAccess.get().set(1);

					System.out.println(String.format(
							"%d, %d, %d",
							ribosomesVolumeAccess.getLongPosition(0),
							ribosomesVolumeAccess.getLongPosition(1),
							ribosomesVolumeAccess.getLongPosition(2)));
				});

		N5Utils.save(ribosomesVolume, n5, "/volumes/labels/ribosomes", new int[] {64, 64, 64}, new GzipCompression());
		n5.setAttribute("/volumes/labels/ribosomes", "offset", fLabelsOffset);
		n5.setAttribute("/volumes/labels/ribosomes", "resolution", labelsResolution);
	}
}

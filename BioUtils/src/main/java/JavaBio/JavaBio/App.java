package JavaBio.JavaBio;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.mmtf.MmtfActions;
import org.biojava.nbio.structure.io.mmtf.MmtfUtils;

/**
 *
 */
public class App 
{
	public static void main( String[] args ) throws IOException, StructureException
	{
		if ( args.length != 2 ){
			printHelpAndExit();
		}
		File f = new File(args[0]);
		if ( ! f.exists() ) {
			System.err.println("Unable to find file: " + f.getAbsolutePath());
			printHelpAndExit();
		}
		String outputFmt = args[1];
		if ( !validOutputFormat(outputFmt) ){
			System.err.println("Invalid output format: " + outputFmt);
			printHelpAndExit();
		}

		AtomCache cache = new AtomCache();

		Structure s = StructureIO.getStructure(f.toURI().toURL().toString());

		File output = new File(f.getParent(), f.getName() + "." + outputFmt);
		String out = null;

		if ( outputFmt.equalsIgnoreCase("mmtf")){
			MmtfUtils.setUpBioJava();
			MmtfActions.writeToFile(s, output.toPath());
		}else{
			if ( outputFmt.equalsIgnoreCase("cif")){
				out = s.toMMCIF();
			}else{
				out = s.toPDB();
			}
			Files.write(output.toPath(), out.getBytes());
		}

	}

	private static boolean validOutputFormat(String outputFmt) {
		List<String> validTypes = Arrays.asList(new String[]{"pdb", "cif", "mmtf"});
		return validTypes.contains(outputFmt.toLowerCase());
	}

	private static void printHelpAndExit() {
		System.err.println("Usage: java -jar BioUtils.jar $FILE $TYPE, where $FILE is a PDB or mmCIF file you wish to convert and " +
				"$TYPE is the type of the output file [PDB, CIF, MMTF].");
		System.exit(1);
	}
}

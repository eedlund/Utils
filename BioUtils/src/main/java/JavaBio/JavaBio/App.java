package JavaBio.JavaBio;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.StructureIO.StructureFiletype;
import org.biojava.nbio.structure.align.util.AtomCache;

/**
 *
 */
public class App 
{
    public static void main( String[] args ) throws IOException, StructureException
    {
    	if ( args.length != 1 ){
    		printHelpAndExit();
    	}
        File f = new File(args[0]);
        if ( ! f.exists() ) {
        	System.err.println("Unable to find file: " + f.getAbsolutePath());
        	printHelpAndExit();
        }
        StructureFiletype type = StructureIO.guessFiletype(f.getAbsolutePath());
        
        AtomCache cache = new AtomCache();
        
        cache.setUseMmCif(type == StructureFiletype.PDB);
        Structure s = StructureIO.getStructure(f.toURI().toURL().toString());
        
        String out = null;
        String ext = "";
        if ( type == StructureFiletype.PDB){
        	out = s.toMMCIF();
        	ext = ".cif";
        }else{
        	out = s.toPDB();
        	ext = ".pdb";
        }
        File output = new File(f.getParent(), f.getName() + ext);
        Files.write(output.toPath(), out.getBytes());
    }

	private static void printHelpAndExit() {
		System.err.println("Usage: java BioUtils.jar $FILE, where $FILE is a PDB or mmCIF file you wish to convert.");
		System.exit(1);
	}
}

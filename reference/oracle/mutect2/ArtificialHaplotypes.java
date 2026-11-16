// Test-only projection in BAM file order; include every stage-determined tag.
import htsjdk.samtools.*;
import java.io.*;
public class ArtificialHaplotypes {
  public static void main(String[] args) throws Exception {
    try (SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(args[0]));
         PrintWriter out = new PrintWriter(new BufferedWriter(new OutputStreamWriter(System.out)))) {
      for (SAMRecord r : reader) {
        if (!"ArtificialHaplotypeRG".equals(r.getStringAttribute("RG"))) continue;
        out.println(r.getReferenceName()+"\t"+r.getAlignmentStart()+"\t"+r.getCigarString()+"\t"+r.getReadString()
          +"\tHC:i:"+r.getIntegerAttribute("HC")+"\tRG:Z:"+r.getStringAttribute("RG")+"\tCR:Z:"+r.getStringAttribute("CR"));
        for (var tag : r.getAttributes()) if (!java.util.Set.of("HC","RG","CR").contains(tag.tag))
          throw new IllegalStateException("Unexpected artificial tag: "+tag.tag);
        if (!r.getReadName().startsWith("HC_") || r.getFlags()!=16) throw new IllegalStateException("Unexpected artificial record: "+r);
      }
    }
  }
}

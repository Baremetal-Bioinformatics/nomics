// Test oracle only: full stage-three biological record projection in file order.
// Includes AR/CR and genotype-dependent XA, preserving pre-existing XA text.
import htsjdk.samtools.*;
import java.io.*;
import java.util.*;
public class BiologicalReads {
  static String value(Object v) {
    if (v instanceof Number n) {
      if (v instanceof Float f) return "f:"+String.format("%08x",Float.floatToRawIntBits(f));
      return "i:"+n.longValue();
    }
    if (v instanceof String s) return "Z:"+s;
    if (v instanceof Character c) return "A:"+c;
    throw new IllegalStateException("Unsupported tag value type: "+v.getClass());
  }
  public static void main(String[] args) throws Exception {
    try (SamReader reader=SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(args[0]));
         PrintWriter out=new PrintWriter(new BufferedWriter(new OutputStreamWriter(System.out)))) {
      for (SAMRecord r:reader) {
        if ("ArtificialHaplotypeRG".equals(r.getStringAttribute("RG"))) continue;
        out.print(r.getReadName()+"\t"+r.getFlags()+"\t"+r.getReferenceName()+"\t"+r.getAlignmentStart()+"\t"+r.getMappingQuality()+"\t"+r.getCigarString()+"\t"+r.getMateReferenceName()+"\t"+r.getMateAlignmentStart()+"\t"+r.getInferredInsertSize()+"\t"+r.getReadString()+"\t"+HexFormat.of().formatHex(r.getBaseQualities()));
        var tags=new ArrayList<>(r.getAttributes()); tags.sort(Comparator.comparing(t->t.tag));
        for(var t:tags) out.print("\t"+t.tag+":"+value(t.value));
        out.println();
      }
    }
  }
}

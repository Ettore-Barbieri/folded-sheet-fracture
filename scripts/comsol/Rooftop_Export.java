/*
 * Rooftop_Export.java  --  export the displacement field of saved models,
 * for the 3D deformation figures. Loads the .mph files written by the
 * Rooftop_FE_* runs (each holds its LAST case) and writes, for each,
 * a text file with x, y, z, u, v, w at the mesh nodes (state 1: tip at c2),
 * plus a small parameter file.
 *
 * Run from scripts/comsol:
 *   /Applications/COMSOL64/Multiphysics/bin/comsol compile Rooftop_Export.java
 *   /Applications/COMSOL64/Multiphysics/bin/comsol batch -inputfile Rooftop_Export.class
 */
import com.comsol.model.*;
import com.comsol.model.util.*;
import java.io.*;

public class Rooftop_Export {
  // Absolute paths: COMSOL's security settings block reading the working directory.
  static final String[] DIRS = {
    "../../results/comsol",
    "/path/to/folded-sheet-fracture/results/comsol",
    "/path/to/folded-sheet-fracture/results/comsol"
  };
  static final String[] FILES = {"Rooftop_FE_ModeI_a10", "Rooftop_FE_ModeII", "Rooftop_FE_ModeIII"};
  static final String[] PARS = {"alpha0", "alpha_g", "lc", "a_over_lc", "a_crack", "delta_over_lc", "delta_tip", "W0", "L0"};

  public static void main(String[] args) { run(); }

  public static Model run() {
    Model last = null;
    String DIR = null;
    for (String d : DIRS) { try { if (new File(d).isDirectory()) { DIR = d; break; } } catch (Exception e) { } }
    if (DIR == null) { System.out.println("No results folder found."); return null; }
    for (String f : FILES) {
      try {
        String path = DIR + File.separator + f + ".mph";
        System.out.println("Loading " + path);
        Model m = ModelUtil.load("M_" + f, path);
        String out = DIR + File.separator + "deform_" + f;
        try { m.result().export().remove("expDef"); } catch (Exception e) { }
        m.result().export().create("expDef", "Data");
        m.result().export("expDef").set("data", "dset1");
        m.result().export("expDef").set("expr", new String[]{"u", "v", "w"});
        m.result().export("expDef").set("filename", out + ".txt");
        m.result().export("expDef").run();
        PrintWriter w = new PrintWriter(new FileWriter(out + "_params.txt"));
        for (String p : PARS) {
          try { w.println(p + " " + m.param().evaluate(p)); } catch (Exception e) { w.println(p + " NaN"); }
        }
        w.close();
        System.out.println("Wrote " + out + ".txt and _params.txt");
        last = m;
      } catch (Exception e) {
        System.out.println("FAILED for " + f + ": " + e.getMessage());
      }
    }
    return last;
  }
}

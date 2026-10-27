/*
 * Rooftop_FE_ModeI_Slope.java  --  COMSOL 6.4 check of the rooftop crack paper.
 * THIS FILE: Mode I, a/lc = 10: facet slope alpha = 0.02 ... 1 at fixed lc/h (limits of the shallow-shell theory).
 * (All Mode I files share the same code; only the switches below differ.)
 *
 * Problem: a symmetric rooftop fold (ridge on x = 0, facets z = -alpha|x|),
 * a crack on y = 0 from x = c1 = -delta-2a to x = c2 = -delta, and a
 * remote membrane tension N_inf along y (parallel to the ridge).
 *
 * Model: half of the sheet (y >= 0), Shell interface, linear, sharp fold.
 *   - Mirror symmetry on the uncracked parts of y = 0:
 *       v = 0   (no y-displacement)
 *       ary = 0 (the shell normal stays in the x-z plane)
 *     "ar" is COMSOL's shell-normal displacement field (arx, ary, arz).
 *   - Crack faces on y = 0 are free.
 *   - Edge load N_inf (N/m) in +y on the edges at y = L0.
 *   - Rigid motion removed: u = w = 0 at the ridge point (0,0,0),
 *     and w = 0 at the far point (W0, 0, -alpha W0). The load is
 *     self-balanced, so these supports carry no force.
 *
 * How G is computed. A short segment of length dx just ahead of the near
 * tip (from c2 to c2+dx) carries its own symmetry constraint.
 *   State 1: segment constrained (tip at c2).
 *   State 2: segment released    (tip at c2+dx).
 * Both states use the SAME mesh, so the difference is clean.
 *   Method 1 (energy):  G = 2 (U2 - U1) / (h dx),
 *       U = (1/2) * integral of N_inf * v along the loaded edge (half model).
 *   Method 2 (closure / VCCT): G = -(1/(h dx)) * sum over segment nodes of
 *       reacf_1 * displacement_2, split into
 *       Gm (from v, membrane) and Gb (from ary, bending).
 *   For a linear model, Method 2 total equals Method 1 exactly, so the two
 *   agreeing checks the post-processing; the value of Method 2 is the split.
 *
 * Calibration: for each alpha, a flat sheet (geometry slope 0) is run with
 * the same mesh rules. G_FE_flat / G_flat shows the finite-domain + mesh error.
 *
 * Output: ModeI_slope.csv and Rooftop_FE_ModeI_Slope.mph in results/comsol (or the folder you run from).
 *
 * Run (macOS, adjust the COMSOL path if needed), from this folder:
 *   /Applications/COMSOL64/Multiphysics/bin/comsol compile Rooftop_FE_ModeI_Slope.java
 *   /Applications/COMSOL64/Multiphysics/bin/comsol batch -inputfile Rooftop_FE_ModeI_Slope.class
 */

import com.comsol.model.*;
import com.comsol.model.util.*;
import java.io.*;
import java.util.*;

public class Rooftop_FE_ModeI_Slope {

  // ---------------- case list ----------------
  static final double[] ALPHAS = {0.02, 0.2};
  static double A_OVER_LC = 10.0;          // changed per case in A100_TEST
  static final double[] DELTAS = {30, 10, 5, 2.28, 1, 0.3, 0.1, 0.03};   // delta / lc
  static double DOMAIN = 25.0;                // W0 = L0 = DOMAIN * a (changed in REFINE_TEST)
  // true: run one case (alpha = 0.02, delta = lc) and hand the model to the GUI
  // for inspection. Set to false for the full sweep once one case works.
  static final boolean ONE_CASE = false;
  // Shear test: alpha = 0.2 only, shear factors 1 and 100, a few deltas.
  static final boolean SHEAR_TEST = false;
  // a/lc test: alpha = 0.02, a/lc = 100, full delta range (theory curve in the paper's G(delta) figure).
  static final boolean A100_TEST = false;
  // Decoupling test: alpha = 0.2, bending (and shear) stiffness x1 and x100.
  // x100 makes lc ten times larger at the same slope, so h/lc ~ 0.033 as for
  // alpha = 0.02. Separates "slope" (geometry) from "thickness" (h/lc) effects.
  static final boolean DECOUPLE_TEST = true;   // used here for the slope test (needs the Section Stiffness feature)
  static final double[] SLOPES = {0.02, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0};
  static final double[] SLOPE_DELTAS = {30, 5, 2.28, 1, 0.3, 0.1, 0.03};
  // Refinement test: alpha = 0.02, a/lc = 10. Is the ~1.3% FE-theory gap from
  // the mesh, the size of the released segment, or the size of the domain?
  static final boolean REFINE_TEST = false;
  static double meshFac = 1;   // multiplies every element size
  static double segFac = 1;    // multiplies the released segment length dx
  static final String CSV_NAME = "ModeI_slope.csv";
  static double bendFac = 1;
  static final double[] A100_DELTAS_FULL = {300, 100, 30, 12.8, 5, 2, 0.6, 0.2, 0.06, 0.03};
  // Quick test first: flat case + one delta, to see mesh size and time per solve.
  static final boolean A100_QUICK = false;
  static final double[] A100_DELTAS = A100_QUICK ? new double[]{1.0} : A100_DELTAS_FULL;
  static final double[] SHEAR_FACS = {1, 100};
  static final double[] SHEAR_DELTAS = {30, 2.28, 0.3, 0.03};
  static double shearFac = 1;

  // ---------------- material (SI) ----------------
  static final double E = 200e9, NU = 0.3, H = 1e-3, NINF = 1e4;

  static Model model;
  static boolean solsExist = false;

  public static void main(String[] args) { run(); }

  public static Model run() {
    model = ModelUtil.create("Model");
    String outDir = findOutDir();          // null if file access is blocked
    if (outDir != null) {
      try { model.modelPath(outDir); } catch (Exception e) { }
      say("Output folder: " + outDir);
    } else {
      say("No output folder available (COMSOL file access may be restricted). "
        + "Results go to the table 'FE results' under Results > Tables.");
    }

    List<String> rows = new ArrayList<String>();
    rows.add("alpha,lc_m,a_over_lc,delta_over_lc,dx_over_lc,deltaMid_over_lc,flat,"
           + "U1,U2,G_energy_over_Gflat,Gm_closure_over_Gflat,Gb_closure_over_Gflat,"
           + "G_closure_over_Gflat,closure_minus_energy_rel,shear_fac,mesh_fac,seg_fac,domain");

    try {
      buildParameters();
      buildGeometry();
      buildSelections();
      buildPhysics();
      buildMesh();
      buildStudiesAndCouplings();

      if (REFINE_TEST) {
        A_OVER_LC = 10.0;
        // {mesh factor, segment factor, domain}
        double[][] variants = {{1, 1, 25}, {0.5, 1, 25}, {0.25, 1, 25}, {1, 0.5, 25}, {1, 1, 50}, {1, 1, 100}, {1, 1, 200}};
        for (double[] v : variants) {
          meshFac = v[0]; segFac = v[1]; DOMAIN = v[2];
          model.param().set("seg_fac", String.valueOf(segFac));
          model.param().set("W0", DOMAIN + "*a_crack");
          model.param().set("L0", DOMAIN + "*a_crack");
          say("REFINE_TEST: mesh x" + meshFac + ", segment x" + segFac + ", domain " + DOMAIN + "a");
          rows.add(runCase(0.02, 1.0, true));
          writeCsv(outDir, rows);
          for (double d : new double[]{30, 2.28, 1, 0.03}) {
            rows.add(runCase(0.02, d, false));
            writeCsv(outDir, rows);
          }
        }
      } else if (DECOUPLE_TEST) {
        // SLOPE TEST (limits of the shallow-shell theory). The facet slope alpha
        // grows from 0.02 to 1 (fold half-angle up to 45 deg). Bending and shear
        // stiffness are multiplied by (alpha/0.02)^2, so lc stays at its alpha = 0.02
        // value (lc/h ~ 30) and the ratio of shear to bending stiffness, which sets
        // the Kirchhoff limit, is unchanged. Only the slope changes.
        A_OVER_LC = 10.0;
        for (double al : SLOPES) {
          double f = (al / 0.02) * (al / 0.02);
          bendFac = f; shearFac = f;
          model.param().set("bend_fac", String.valueOf(f));
          model.param().set("shear_fac", String.valueOf(f));
          say("SLOPE_TEST: alpha = " + al + ", bending and shear stiffness x" + f);
          rows.add(runCase(al, 1.0, true));
          writeCsv(outDir, rows);
          for (double d : SLOPE_DELTAS) {
            rows.add(runCase(al, d, false));
            writeCsv(outDir, rows);
          }
        }
      } else if (A100_TEST) {
        // a/lc = 100 (full curve), a/lc = 10 (cross-check on the new mesh),
        // a/lc = 1000 (a few points, for the shielding law).
        double[] aList = A100_QUICK ? new double[]{100} : new double[]{100, 10, 1000};
        for (double aL : aList) {
          A_OVER_LC = aL;
          double[] ds = (aL == 100) ? A100_DELTAS
                      : (aL == 10)  ? new double[]{1.0, 0.03}
                      :               new double[]{1000, 134, 10, 1.0, 0.03};
          say("A100_TEST mode: alpha = 0.02, a/lc = " + aL);
          rows.add(runCase(0.02, 1.0, true));
          writeCsv(outDir, rows);
          for (double d : ds) {
            rows.add(runCase(0.02, d, false));
            writeCsv(outDir, rows);
          }
        }
      } else if (SHEAR_TEST) {
        say("SHEAR_TEST mode: alpha = 0.2, shear factors 1 and 100.");
        for (double f : SHEAR_FACS) {
          shearFac = f;
          model.param().set("shear_fac", String.valueOf(f));
          rows.add(runCase(0.2, 1.0, true));
          writeCsv(outDir, rows);
          for (double d : SHEAR_DELTAS) {
            rows.add(runCase(0.2, d, false));
            writeCsv(outDir, rows);
          }
        }
      } else if (ONE_CASE) {
        say("ONE_CASE mode: running a single case. Set ONE_CASE = false for the full sweep.");
        rows.add(runCase(0.02, 1.0, false));
        writeCsv(outDir, rows);
      } else {
        for (double alpha : ALPHAS) {
          // flat-sheet calibration (geometry slope 0, same alpha for lc)
          rows.add(runCase(alpha, 1.0, true));
          writeCsv(outDir, rows);
          for (double d : DELTAS) {
            rows.add(runCase(alpha, d, false));
            writeCsv(outDir, rows);
          }
        }
      }
    } catch (Exception e) {
      // Do not throw: the GUI would discard the model and the log with it.
      say("STOPPED WITH ERROR: " + firstLine(e));
      String m = String.valueOf(e.getMessage());
      String[] lines = m.split("\n");
      for (int i = 1; i < Math.min(lines.length, 12); i++) say("   | " + lines[i].trim());
    }

    if (outDir != null) try {
      model.save(outDir + File.separator + "Rooftop_FE_ModeI_Slope.mph");
      say("Saved Rooftop_FE_ModeI_Slope.mph (last case).");
    } catch (Exception e) { say("Could not save .mph: " + firstLine(e)); }
    flushLog();
    return model;
  }

  // =====================================================================
  static void buildParameters() {
    model.param().set("E0", "200[GPa]", "Young's modulus");
    model.param().set("nu0", "0.3", "Poisson's ratio");
    model.param().set("h0", "1[mm]", "Thickness h");
    model.param().set("alpha0", "0.02", "Facet slope alpha (sets lc)");
    model.param().set("alpha_g", "alpha0", "Slope used in the geometry (0 = flat calibration)");
    model.param().set("bend_fac", "1", "Multiplier on bending stiffness (Section Stiffness only)");
    model.param().set("lc", "sqrt(bend_fac)*h0/(alpha0*sqrt(3*(1-nu0^2)))", "Fold length scale lc");
    model.param().set("a_over_lc", String.valueOf(A_OVER_LC), "a / lc");
    model.param().set("a_crack", "a_over_lc*lc", "Crack half-length a");
    model.param().set("delta_over_lc", "1", "delta / lc");
    model.param().set("delta_tip", "delta_over_lc*lc", "Near tip to ridge distance");
    model.param().set("seg_fac", "1", "Multiplier on the released segment length");
    model.param().set("dx_seg", "seg_fac*min(lc/20, delta_tip/5)", "Released segment length");
    model.param().set("W0", DOMAIN + "*a_crack", "Facet width (each side)");
    model.param().set("L0", DOMAIN + "*a_crack", "Half-model height");
    model.param().set("ridge_len", "min(5*a_crack, 50*lc)", "Refined length of the ridge (from y = 0)");
    model.param().set("N_inf", "1e4[N/m]", "Remote membrane tension");
    model.param().set("G_flat", "pi*N_inf^2*a_crack/(E0*h0^2)", "Griffith G, flat sheet");
    model.param().set("shear_fac", "1", "Multiplier on transverse shear stiffness (1 = real shell)");
  }

  // =====================================================================
  static void buildGeometry() {
    model.component().create("comp1", true);
    model.component("comp1").geom().create("geom1", 3);

    // Profile on y = 0, drawn in the xz work plane (local x = x, local y = z).
    model.component("comp1").geom("geom1").create("wp1", "WorkPlane");
    model.component("comp1").geom("geom1").feature("wp1").set("quickplane", "xz");
    model.component("comp1").geom("geom1").feature("wp1").geom().create("pol1", "Polygon");
    model.component("comp1").geom("geom1").feature("wp1").geom().feature("pol1").set("source", "table");
    model.component("comp1").geom("geom1").feature("wp1").geom().feature("pol1").set("type", "open");
    model.component("comp1").geom("geom1").feature("wp1").geom().feature("pol1").set("table", new String[][]{
      {"-W0",                    "-alpha_g*W0"},
      {"-delta_tip-2*a_crack",   "-alpha_g*(delta_tip+2*a_crack)"},
      {"0",                      "0"},
      {"W0",                     "-alpha_g*W0"}
    });

    // The xz work plane has its normal along -y, so a plain extrude goes to
    // y < 0. "reverse" sends it to y in [0, L0].
    model.component("comp1").geom("geom1").create("ext1", "Extrude");
    model.component("comp1").geom("geom1").feature("ext1").set("workplane", "wp1");
    model.component("comp1").geom("geom1").feature("ext1").selection("input").set(new String[]{"wp1"});
    // Two layers: the ridge is split at y = ridge_len, so only its part near
    // the crack gets the fine ridge mesh.
    model.component("comp1").geom("geom1").feature("ext1").setIndex("distance", "ridge_len", 0);
    model.component("comp1").geom("geom1").feature("ext1").setIndex("distance", "L0", 1);
    model.component("comp1").geom("geom1").feature("ext1").set("reverse", true);
    // The two tip points (c2 and c2+dx) are added as single points on the
    // y = 0 edge. Putting them in the extruded profile would create strips
    // dx wide running the full height of the sheet, which explodes the mesh.
    String[][] pts = {{"pt1", "-delta_tip", "-alpha_g*delta_tip"},
                      {"pt2", "-delta_tip+dx_seg", "-alpha_g*(delta_tip-dx_seg)"}};
    for (String[] q : pts) {
      model.component("comp1").geom("geom1").create(q[0], "Point");
      model.component("comp1").geom("geom1").feature(q[0]).setIndex("p", q[1], 0);
      model.component("comp1").geom("geom1").feature(q[0]).setIndex("p", "0", 1);
      model.component("comp1").geom("geom1").feature(q[0]).setIndex("p", q[2], 2);
    }
    model.component("comp1").geom("geom1").run();
  }

  // =====================================================================
  // Box selections. Bounds are numbers, refreshed for every case.
  static final String[] BOX1 = {"selFar", "selCrack", "selSeg", "selNear", "selRight", "selTop", "selRidge"};
  static final String[] BOX0 = {"selTipPts", "selOrigin", "selFarPt"};
  static final String[] BOX2 = {"selFaces"};

  static void buildSelections() {
    for (String t : BOX1) { model.component("comp1").selection().create(t, "Box"); model.component("comp1").selection(t).set("entitydim", 1); model.component("comp1").selection(t).set("condition", "inside"); }
    for (String t : BOX0) { model.component("comp1").selection().create(t, "Box"); model.component("comp1").selection(t).set("entitydim", 0); model.component("comp1").selection(t).set("condition", "inside"); }
    for (String t : BOX2) { model.component("comp1").selection().create(t, "Box"); model.component("comp1").selection(t).set("entitydim", 2); model.component("comp1").selection(t).set("condition", "inside"); }
    model.component("comp1").selection().create("selSym", "Union");
    model.component("comp1").selection("selSym").set("entitydim", 1);
    model.component("comp1").selection("selSym").set("input", new String[]{"selFar", "selNear", "selRight"});
  }

  static void setBox(String tag, double x0, double x1, double y0, double y1) {
    model.component("comp1").selection(tag).set("xmin", String.valueOf(x0));
    model.component("comp1").selection(tag).set("xmax", String.valueOf(x1));
    model.component("comp1").selection(tag).set("ymin", String.valueOf(y0));
    model.component("comp1").selection(tag).set("ymax", String.valueOf(y1));
    model.component("comp1").selection(tag).set("zmin", "-1e4");
    model.component("comp1").selection(tag).set("zmax", "1e4");
  }

  static void updateSelections(double lc, double a, double delta, double dx) {
    double W0 = DOMAIN * a, L0 = DOMAIN * a;
    double c1 = -delta - 2 * a, c2 = -delta;
    double t = 0.01 * dx;               // tolerance, far below the shortest edge
    double tL = 1e-6 * L0;
    setBox("selFar",   -W0 - t, c1 + t,       -t, t);
    setBox("selCrack",  c1 - t, c2 + t,       -t, t);
    setBox("selSeg",    c2 - t, c2 + dx + t,  -t, t);
    setBox("selNear",   c2 + dx - t, t,       -t, t);
    setBox("selRight", -t, W0 + t,            -t, t);
    setBox("selTop",   -W0 - tL, W0 + tL,     L0 - tL, L0 + tL);
    setBox("selRidge", -t, t,                 -t, Math.min(5 * a, 50 * lc) + tL);   // near part of the ridge only
    setBox("selTipPts", c2 - t, c2 + dx + t,  -t, t);
    setBox("selOrigin", -t, t,                -t, t);
    setBox("selFarPt",  W0 - tL, W0 + tL,     -t, t);
    setBox("selFaces", -W0 - tL, W0 + tL,     -tL, L0 + tL);   // all five facet strips
  }

  // =====================================================================
  static void buildPhysics() {
    model.component("comp1").material().create("mat1", "Common");
    model.component("comp1").material("mat1").propertyGroup("def").set("youngsmodulus", "E0");
    model.component("comp1").material("mat1").propertyGroup("def").set("poissonsratio", "nu0");
    model.component("comp1").material("mat1").propertyGroup("def").set("density", "1000[kg/m^3]");

    model.component("comp1").physics().create("shell", "Shell", "geom1");
    model.component("comp1").physics("shell").feature("to1").set("d", "h0");

    // Remote tension
    model.component("comp1").physics("shell").create("el1", "EdgeLoad", 1);
    model.component("comp1").physics("shell").feature("el1").selection().named("selTop");
    model.component("comp1").physics("shell").feature("el1").set("forceReferenceLength", new String[]{"0", "N_inf", "0"});

    // Symmetry on the ligaments, and on the short segment (switched off in state 2)
    addSymmetry("pcV", "pcA", "selSym");
    addSymmetry("pcSegV", "pcSegA", "selSeg");

    // Section Stiffness: same membrane and bending stiffness as the material,
    // transverse shear stiffness set by hand to shear_fac * kappa*G*h, kappa = 5/6.
    // shear_fac = 1 reproduces the ordinary shell; shear_fac = 100 is close to Kirchhoff.
    model.component("comp1").physics("shell").create("ses1", "SectionStiffness", 2);
    model.component("comp1").physics("shell").feature("ses1").selection().all();
    if (!SHEAR_TEST && !DECOUPLE_TEST) model.component("comp1").physics("shell").feature("ses1").active(false);
    // "From material" would need these matrices stored in the material, so all
    // are given by hand. Isotropic plate, Voigt order (11, 22, 12), full 3x3.
    String A = "E0*h0/(1-nu0^2)", Dm = "bend_fac*E0*h0^3/(12*(1-nu0^2))";
    model.component("comp1").physics("shell").feature("ses1").set("DA_mat", "userdef");
    model.component("comp1").physics("shell").feature("ses1").set("DA", iso3(A));
    model.component("comp1").physics("shell").feature("ses1").set("DD_mat", "userdef");
    model.component("comp1").physics("shell").feature("ses1").set("DD", iso3(Dm));
    model.component("comp1").physics("shell").feature("ses1").set("DB_mat", "userdef");
    model.component("comp1").physics("shell").feature("ses1").set("DB", new String[]{"0","0","0","0","0","0","0","0","0"});
    model.component("comp1").physics("shell").feature("ses1").set("DAs_mat", "userdef");
    model.component("comp1").physics("shell").feature("ses1").set("DAs", new String[]{
      "shear_fac*(5/6)*E0/(2*(1+nu0))*h0", "0", "0", "shear_fac*(5/6)*E0/(2*(1+nu0))*h0"});

    // Rigid-body supports
    model.component("comp1").physics("shell").create("disp1", "Displacement0", 0);
    model.component("comp1").physics("shell").feature("disp1").selection().named("selOrigin");
    model.component("comp1").physics("shell").feature("disp1").set("Direction", new String[]{"prescribed", "free", "prescribed"});
    model.component("comp1").physics("shell").create("disp2", "Displacement0", 0);
    model.component("comp1").physics("shell").feature("disp2").selection().named("selFarPt");
    model.component("comp1").physics("shell").feature("disp2").set("Direction", new String[]{"free", "free", "prescribed"});
  }

  static String[] iso3(String k) {
    String d = k, o = "nu0*(" + k + ")", g = "(1-nu0)/2*(" + k + ")";
    return new String[]{d, o, "0", o, d, "0", "0", "0", g};
  }

  // Mirror symmetry about y = 0 for a shell: v = 0 and ary = 0.
  // First choice: two pointwise constraints (exact and transparent).
  // Fallback: COMSOL's own shell symmetry feature, under the names it may use.
  // In the fallback, tagA is not created; the toggling code skips it.
  static String symMode = null;
  static final String[] SYM_TYPES = {"Symmetry", "SymmetrySolid", "SymmetryShell", "SymmetryAntisymmetry"};

  static void addSymmetry(String tagV, String tagA, String sel) {
    if (symMode == null || symMode.equals("pointwise")) {
      try {
        model.component("comp1").physics("shell").create(tagV, "PointwiseConstraint", 1);
        model.component("comp1").physics("shell").feature(tagV).selection().named(sel);
        model.component("comp1").physics("shell").feature(tagV).set("constraintExpression", "v");
        model.component("comp1").physics("shell").create(tagA, "PointwiseConstraint", 1);
        model.component("comp1").physics("shell").feature(tagA).selection().named(sel);
        model.component("comp1").physics("shell").feature(tagA).set("constraintExpression", "ary");
        symMode = "pointwise";
        say("Symmetry via pointwise constraints v = 0, ary = 0 on " + sel);
        return;
      } catch (Exception e) {
        say("PointwiseConstraint not accepted (" + e.getMessage() + "); trying shell Symmetry feature.");
        try { model.component("comp1").physics("shell").feature().remove(tagV); } catch (Exception e2) { }
        try { model.component("comp1").physics("shell").feature().remove(tagA); } catch (Exception e2) { }
      }
    }
    for (String type : SYM_TYPES) {
      try {
        model.component("comp1").physics("shell").create(tagV, type, 1);
        model.component("comp1").physics("shell").feature(tagV).selection().named(sel);
        symMode = type;
        say("Symmetry via shell feature type '" + type + "' on " + sel);
        return;
      } catch (Exception e) {
        try { model.component("comp1").physics("shell").feature().remove(tagV); } catch (Exception e2) { }
      }
    }
    throw new RuntimeException("No symmetry feature could be created. Add one Symmetry edge "
      + "condition in the Shell GUI, save as Java, and send me the line that creates it.");
  }

  static void setSegActive(boolean on) {
    model.component("comp1").physics("shell").feature("pcSegV").active(on);
    if ("pointwise".equals(symMode)) model.component("comp1").physics("shell").feature("pcSegA").active(on);
  }

  // =====================================================================
  static void buildMesh() {
    model.component("comp1").mesh().create("mesh1");
    model.component("comp1").mesh("mesh1").feature("size").set("custom", "on");
    String[][] local = {
      {"szTip",   "0", "selTipPts"},
      {"szSeg",   "1", "selSeg"},
      {"szNear",  "1", "selNear"},
      {"szCrack", "1", "selCrack"},
      {"szRidge", "1", "selRidge"}
    };
    for (String[] s : local) {
      model.component("comp1").mesh("mesh1").create(s[0], "Size");
      model.component("comp1").mesh("mesh1").feature(s[0]).selection().geom("geom1", Integer.parseInt(s[1]));
      model.component("comp1").mesh("mesh1").feature(s[0]).selection().named(s[2]);
      model.component("comp1").mesh("mesh1").feature(s[0]).set("custom", "on");
      model.component("comp1").mesh("mesh1").feature(s[0]).set("hmaxactive", true);
      model.component("comp1").mesh("mesh1").feature(s[0]).set("hgradactive", true);
      model.component("comp1").mesh("mesh1").feature(s[0]).set("hgrad", "1.15");
    }
    // The surface mesher needs an explicit face selection when created from
    // code; left empty it meshes nothing (this was the "not meshed" error).
    model.component("comp1").mesh("mesh1").create("ftri1", "FreeTri");
    model.component("comp1").mesh("mesh1").feature("ftri1").selection().geom("geom1", 2);
    model.component("comp1").mesh("mesh1").feature("ftri1").selection().named("selFaces");
  }

  static void updateMesh(double lc, double a, double delta, double dx) {
    double f = meshFac;
    model.component("comp1").mesh("mesh1").feature("size").set("hmax", String.valueOf(f * a));
    model.component("comp1").mesh("mesh1").feature("size").set("hmin", String.valueOf(f * dx / 20));
    model.component("comp1").mesh("mesh1").feature("size").set("hgrad", "1.2");
    model.component("comp1").mesh("mesh1").feature("szTip").set("hmax",   String.valueOf(f * dx / 4));
    model.component("comp1").mesh("mesh1").feature("szSeg").set("hmax",   String.valueOf(f * dx / 4));
    model.component("comp1").mesh("mesh1").feature("szNear").set("hmax",  String.valueOf(f * Math.min(lc / 4, Math.max(dx / 4, (delta - dx) / 4))));
    model.component("comp1").mesh("mesh1").feature("szCrack").set("hmax", String.valueOf(f * a / 20));
    model.component("comp1").mesh("mesh1").feature("szRidge").set("hmax", String.valueOf(f * lc / 2));
    model.component("comp1").mesh("mesh1").run();
  }

  // =====================================================================
  static void buildStudiesAndCouplings() {
    // U (half model) = 1/2 * integral of N_inf * v over the loaded edges
    model.component("comp1").cpl().create("intTop", "Integration");
    model.component("comp1").cpl("intTop").selection().geom("geom1", 1);
    model.component("comp1").cpl("intTop").selection().named("selTop");
    model.component("comp1").cpl("intTop").set("opname", "intTop");
    // Sum over mesh nodes of the segment (for reaction forces)
    model.component("comp1").cpl().create("sumSeg", "Integration");
    model.component("comp1").cpl("sumSeg").selection().geom("geom1", 1);
    model.component("comp1").cpl("sumSeg").selection().named("selSeg");
    model.component("comp1").cpl("sumSeg").set("opname", "sumSeg");
    model.component("comp1").cpl("sumSeg").set("method", "summation");

    model.study().create("std1");
    model.study("std1").create("stat", "Stationary");
    model.study().create("std2");
    model.study("std2").create("stat", "Stationary");
  }

  static void solve(String std, String sol, String dset) {
    if (solsExist) {
      try { model.result().dataset().remove(dset); } catch (Exception e) { }
      try { model.sol().remove(sol); } catch (Exception e) { }
    }
    model.sol().create(sol);
    model.sol(sol).study(std);
    model.sol(sol).attach(std);
    model.sol(sol).createAutoSequence(std);
    forceDirect(sol);
    model.sol(sol).runAll();
    try { model.result().dataset().remove(dset); } catch (Exception e) { }
    model.result().dataset().create(dset, "Solution");
    model.result().dataset(dset).set("solution", sol);
  }

  // For large models COMSOL's automatic setup picks an iterative solver,
  // which did not converge here. Force the direct solver (PARDISO) instead.
  static void forceDirect(String sol) {
    try {
      model.sol(sol).feature("s1").feature("dDef").set("linsolver", "pardiso");
      // At a/lc = 1000 the sheet is ~750 m wide and 1 mm thick, so the matrix
      // is badly conditioned. The residual is ~1e-13, but COMSOL's error
      // estimate (based on the condition number) exceeds its limit and the
      // solve is rejected. Accuracy is checked instead by energy vs closure,
      // the flat calibration and the theory curve.
      try { model.sol(sol).feature("s1").feature("dDef").set("errorchk", "off"); }
      catch (Exception e2) { say("   (could not switch off the error check: " + firstLine(e2) + ")"); }
      model.sol(sol).feature("s1").feature("fc1").set("linsolver", "dDef");
    } catch (Exception e) {
      say("   (could not force the direct solver: " + firstLine(e) + ")");
    }
  }

  static double eval(String dset, String expr) {
    String t = "gevTmp";
    try { model.result().numerical().remove(t); } catch (Exception e) { }
    model.result().numerical().create(t, "EvalGlobal");
    model.result().numerical(t).set("data", dset);
    model.result().numerical(t).set("expr", new String[]{expr});
    double[][] r = model.result().numerical(t).getReal();
    return r[0][0];
  }

  // =====================================================================
  static String runCase(double alpha, double dOverLc, boolean flat) {
    double lc = Math.sqrt(bendFac) * H / (alpha * Math.sqrt(3 * (1 - NU * NU)));
    double a = A_OVER_LC * lc;
    double delta = dOverLc * lc;
    double dx = segFac * Math.min(lc / 20, delta / 5);
    model.param().set("a_over_lc", String.valueOf(A_OVER_LC));
    double Gflat = Math.PI * NINF * NINF * a / (E * H * H);

    say(String.format("---- alpha=%.3g  delta/lc=%.4g  %s  (lc=%.4g m, dx/lc=%.3g)",
        alpha, dOverLc, flat ? "FLAT CALIBRATION" : "", lc, dx / lc));

    model.param().set("alpha0", String.valueOf(alpha));
    model.param().set("alpha_g", flat ? "0" : "alpha0");
    model.param().set("delta_over_lc", String.valueOf(dOverLc));

    model.component("comp1").geom("geom1").run();
    updateSelections(lc, a, delta, dx);
    checkSelections();
    updateMesh(lc, a, delta, dx);
    meshReport();

    // State 2: segment released (tip at c2+dx). Solved first.
    setSegActive(false);
    solve("std2", "sol2", "dset2");
    // State 1: segment constrained (tip at c2). Solved last, so its
    // constraints are the active ones when reaction forces are evaluated.
    setSegActive(true);
    solve("std1", "sol1", "dset1");
    solsExist = true;

    double U1 = energy("dset1");
    double U2 = energy("dset2");
    double Ge = 2 * (U2 - U1) / (H * dx);

    double Gm = Double.NaN, Gb = Double.NaN;
    try {
      // Work of all released reactions: displacements (membrane) and
      // shell-normal displacements (bending). With the pointwise form only
      // v and ary carry reactions; the full sums also cover the fallback.
      Gm = -nodeSum("dset1", "reacf(u)*withsol('sol2',u)+reacf(v)*withsol('sol2',v)+reacf(w)*withsol('sol2',w)") / (H * dx);
      Gb = -nodeSum("dset1", "reacf(arx)*withsol('sol2',arx)+reacf(ary)*withsol('sol2',ary)+reacf(arz)*withsol('sol2',arz)") / (H * dx);
    } catch (Exception e) {
      say("Closure (VCCT) evaluation failed: " + e.getMessage());
    }
    double Gc = Gm + Gb;
    double rel = (Gc - Ge) / Ge;

    if (U1 == U2) say("WARNING: U1 == U2. The segment release did not change the model.");
    say(String.format("   G_energy/G_flat = %.5f   Gm/G_flat = %.5f   Gb/G_flat = %.5f   (closure-energy)/energy = %.2e",
        Ge / Gflat, Gm / Gflat, Gb / Gflat, rel));

    double dMid = (delta - dx / 2) / lc;
    return String.format(Locale.US, "%.4g,%.6e,%.4g,%.6g,%.6g,%.6g,%d,%.10e,%.10e,%.6f,%.6f,%.6f,%.6f,%.3e,%.4g,%.4g,%.4g,%.4g",
        alpha, lc, A_OVER_LC, dOverLc, dx / lc, dMid, flat ? 1 : 0,
        U1, U2, Ge / Gflat, Gm / Gflat, Gb / Gflat, Gc / Gflat, rel, shearFac, meshFac, segFac, DOMAIN);
  }

  // Stop early, with a clear message, if a box selection is empty.
  static void checkSelections() {
    String[] all = {"selFar", "selCrack", "selSeg", "selNear", "selRight", "selTop", "selRidge",
                    "selTipPts", "selOrigin", "selFarPt", "selFaces"};
    int[] expected = {1, 1, 1, 1, 1, 3, 1, 2, 1, 1, 6};
    int[] dims     = {1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 2};
    for (int i = 0; i < all.length; i++) {
      int n;
      try { n = model.component("comp1").selection(all[i]).entities(dims[i]).length; }
      catch (Exception e) { say("(selection check skipped: " + e.getMessage() + ")"); return; }
      if (n != expected[i])
        say("WARNING: selection " + all[i] + " has " + n + " entities, expected " + expected[i]);
      if (n == 0) throw new RuntimeException("Empty selection " + all[i] + " - geometry is not where expected.");
    }
  }

  static void writeCsv(String dir, List<String> rows) {
    updateTable(rows);
    if (dir == null) return;
    try {
      PrintWriter w = new PrintWriter(new FileWriter(dir + File.separator + CSV_NAME));
      for (String r : rows) w.println(r);
      w.close();
    } catch (Exception e) { say("CSV write failed: " + e.getMessage()); }
  }

  // Messages go to the console and to a buffer. At the very end the buffer is
  // written to parameters log001, log002, ... (Global Definitions > Parameters,
  // Description column), so the log is visible in the GUI. Writing them only at
  // the end avoids touching parameters while meshes and solutions are in use.
  static List<String> logBuf = new ArrayList<String>();
  static void say(String s) { System.out.println(s); logBuf.add(s); }

  static void flushLog() {
    for (int i = 0; i < logBuf.size(); i++) {
      try {
        String d = logBuf.get(i).trim().replace('"', '\'');
        if (d.length() > 240) d = d.substring(0, 240);
        model.param().set(String.format("log%03d", i + 1), "0", d);
      } catch (Exception e) { }
    }
  }

  // ---- energy U of the half model, with fallbacks ----
  // 1: coupling operator on the loaded edges, U = 1/2 * int N_inf v dl
  // 2: Line Integration result feature on the same edges
  // 3: total strain energy of the Shell interface (equal to 1 and 2 for this linear problem)
  static int energyRoute = 0;
  static double energy(String dset) {
    String[] names = {"edge operator intTop", "Line Integration", "shell.Ws_tot"};
    for (int r = (energyRoute == 0 ? 1 : energyRoute); r <= 3; r++) {
      try {
        double U;
        if (r == 1) U = eval(dset, "0.5*intTop(N_inf*v)");
        else if (r == 2) U = lineInt(dset, "selTop", "0.5*N_inf*v", false);
        else U = eval(dset, "shell.Ws_tot");
        if (energyRoute != r) say("   energy route used: " + names[r - 1]);
        energyRoute = r;
        return U;
      } catch (Exception e) {
        say("   energy route '" + names[r - 1] + "' failed: " + firstLine(e));
        if (energyRoute == r) energyRoute = 0;
      }
    }
    // Diagnostics: can the dataset evaluate anything at all?
    try { say("   check: global '1' on " + dset + " = " + eval(dset, "1")); }
    catch (Exception e) { say("   check: global '1' on " + dset + " failed: " + firstLine(e)); }
    try { say("   check: max of v on top edges = " + lineInt(dset, "selTop", "v", false)); }
    catch (Exception e) { say("   check: integral of v failed: " + firstLine(e)); }
    throw new RuntimeException("All energy routes failed.");
  }

  // ---- sum over mesh nodes of the released segment, with a fallback ----
  static int sumRoute = 0;
  static double nodeSum(String dset, String expr) {
    for (int r = (sumRoute == 0 ? 1 : sumRoute); r <= 2; r++) {
      try {
        double S = (r == 1) ? eval(dset, "sumSeg(" + expr + ")") : lineInt(dset, "selSeg", expr, true);
        if (sumRoute != r) say("   node-sum route used: " + (r == 1 ? "operator sumSeg" : "Line Integration (summation)"));
        sumRoute = r;
        return S;
      } catch (Exception e) {
        say("   node-sum route " + r + " failed: " + firstLine(e));
        if (sumRoute == r) sumRoute = 0;
      }
    }
    throw new RuntimeException("All node-sum routes failed.");
  }

  static double lineInt(String dset, String sel, String expr, boolean summation) {
    String t = "intTmp";
    try { model.result().numerical().remove(t); } catch (Exception e) { }
    model.result().numerical().create(t, "IntLine");
    model.result().numerical(t).set("data", dset);
    model.result().numerical(t).selection().named(sel);
    if (summation) model.result().numerical(t).set("method", "summation");
    model.result().numerical(t).set("expr", new String[]{expr});
    double[][] r = model.result().numerical(t).getReal();
    return r[0][0];
  }

  static String firstLine(Exception e) {
    String m = String.valueOf(e.getMessage());
    int i = m.indexOf('\n');
    return i > 0 ? m.substring(0, i) : m;
  }

  static void meshReport() {
    int nt = -1;
    try {
      nt = model.component("comp1").mesh("mesh1").getNumElem("tri");
      say("   mesh: " + nt + " triangles, "
        + model.component("comp1").mesh("mesh1").getNumElem("edg") + " edge elements");
    } catch (Exception e) { say("   (mesh report unavailable: " + firstLine(e) + ")"); }
    if (nt == 0) throw new RuntimeException("The mesh is empty: no triangles were generated.");
  }


  // Results table inside the model (Results > Tables > FE results).
  // Always filled, so nothing is lost if file access is blocked.
  // Right-click the table in the GUI to export it.
  static boolean tableMade = false;
  static void updateTable(List<String> rows) {
    try {
      if (!tableMade) {
        model.result().table().create("tblFE", "Table");
        model.result().table("tblFE").label("FE results");
        tableMade = true;
      }
      int nr = rows.size() - 1;
      String[] head = rows.get(0).split(",");
      double[][] data = new double[nr][head.length];
      for (int i = 0; i < nr; i++) {
        String[] f = rows.get(i + 1).split(",");
        for (int j = 0; j < head.length; j++) {
          try { data[i][j] = Double.parseDouble(f[j]); } catch (Exception e) { data[i][j] = Double.NaN; }
        }
      }
      model.result().table("tblFE").setColumnHeaders(head);
      model.result().table("tblFE").setTableData(data);
    } catch (Exception e) { say("Table update failed: " + e.getMessage()); }
  }

  // COMSOL's security settings may block System.getProperty, so the output
  // folder is written here. The first folder that exists is used.
  // Add your own path to the list if you move the project.
  static final String[] OUT_DIRS = {
    "/path/to/folded-sheet-fracture/results/comsol",
    "/path/to/folded-sheet-fracture/results/comsol",
    "../../results/comsol",        // relative to scripts/comsol, where the batch is started
    "."
  };

  static String findOutDir() {
    for (String d : OUT_DIRS) {
      try { if (new File(d).isDirectory()) return d; } catch (Exception e) { }
    }
    return null;
  }
}

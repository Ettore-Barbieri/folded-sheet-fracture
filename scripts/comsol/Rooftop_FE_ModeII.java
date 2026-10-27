/*
 * Rooftop_FE_ModeII.java  --  COMSOL 6.4 check of mode II.
 * (Modes II and III share this code; only MODE differs.) Check of modes II and III
 * (the antisymmetric family) for a crack perpendicular to a rooftop fold.
 *
 * Same rooftop, crack and mesh rules as Rooftop_FE.java (mode I).
 * Ridge on x = 0, facets z = -alpha|x|, crack on y = 0 from
 * x = c1 = -delta-2a to x = c2 = -delta. Half model y >= 0, Shell interface,
 * linear, sharp fold.
 *
 * MODE (fixed in this file to "II"):
 *   "II"  : remote in-plane shear N_xy = S_inf. Crack faces free.
 *           Edge loads: on the top edge y = L0, S_inf along the facet
 *           direction t (+x side); on the side edges x = +-W0, +-S_inf along y.
 *   "III" : no remote load. Both crack faces pushed apart by a uniform
 *           transverse shear V0 (force per unit length) normal to the facet.
 *           In the half model we load the face y = 0+ along +n; antisymmetry
 *           supplies the opposite load on y = 0-. NOTE: here the crack faces
 *           are NOT traction-free (see the mode III note).
 *
 * Antisymmetry on the uncracked parts of y = 0 (u_x and w odd in y, u_y even):
 *     u = 0, w = 0                  (global x and z displacements)
 *     ar . t = 0                    (shell-normal field: its component along the
 *                                    facet direction t = (1,0,s)/sqrt(1+s^2),
 *                                    s = alpha_g*sign(-x), is odd in y)
 *   "ar" is COMSOL's shell-normal displacement field (arx, ary, arz).
 *   In mode I (symmetric) the complementary set v = 0, ary = 0 was used.
 *
 * Rigid motion: the antisymmetry line removes x and z translation and the
 * rotation about y. We add v = 0 at the origin, v = 0 at the far right point
 * on y = 0 and w = 0 at the top of the ridge (0, L0, 0). The loads are
 * self-balanced, so these supports carry no force.
 *
 * G, as in mode I: a segment of length dx ahead of the near tip carries its
 * own antisymmetry constraint. State 1: constrained (tip at c2). State 2:
 * released (tip at c2 + dx). Same mesh in both.
 *   Method 1 (energy):  G = 2 (U2 - U1) / (h dx), U = 1/2 * work of the loads
 *                       on the half model.
 *   Method 2 (closure): G = -(1/(h dx)) * sum over segment nodes of
 *                       reacf_1 . displacement_2, split into
 *                       Gm (in-plane, along t: sliding) and
 *                       Gt (out of plane, along n, plus rotations: tearing).
 *   In mode III the released segment is NOT loaded in state 2. The missing
 *   face-load work on it is of higher order in dx and vanishes as dx -> 0.
 *   Both methods then refer to the same two states and must agree.
 *
 * Calibration: for each case list a flat sheet (geometry slope 0) is run with
 * the same mesh rules; results are divided by it (G/G_flat_FE).
 * Also printed: G_flat_FE / G_flat_theory (mode II: pi S^2 a/(E h^2);
 * mode III: Kirchhoff pi V0^2 a^3/(4 c3 D h), c3 = (1-nu)(3+nu)).
 *
 * Use alpha = 0.02: then the developed lengths differ from the x lengths by
 * a factor 1.0002 and no length conversion is needed at the 0.1% level.
 *
 * Run (macOS), from this folder:
 *   /Applications/COMSOL64/Multiphysics/bin/comsol compile Rooftop_FE_ModeII.java
 *   /Applications/COMSOL64/Multiphysics/bin/comsol batch -inputfile Rooftop_FE_ModeII.class
 * Output: Rooftop_FE_ModeII.csv (or _ModeIII.csv), a log, and the .mph of the
 * last case, in results/comsol.
 */

import com.comsol.model.*;
import com.comsol.model.util.*;
import java.io.*;
import java.util.*;

public class Rooftop_FE_ModeII {

  // ---------------- what to run ----------------
  static final String MODE = "II";
  // QUICK = true: flat + one point (delta = lc) at a/lc = 10, to check set-up and timing.
  static final boolean QUICK = false;
  static final double ALPHA = 0.02;
  static final double[] A_LIST = {10, 100};
  static final double[] DELTAS_A10  = {30, 10, 3, 1, 0.3, 0.1, 0.03};
  static final double[] DELTAS_A100 = {300, 30, 3, 0.3, 0.03};
  static double DOMAIN = 25.0;                     // W0 = L0 = DOMAIN * a
  static double A_OVER_LC = 10.0;

  // ---------------- material and loads (SI) ----------------
  static final double E = 200e9, NU = 0.3, H = 1e-3;
  static final double SINF = 1e4;                  // N/m, mode II
  static final double V0 = 1e2;                    // N/m, mode III

  static Model model;
  static boolean solsExist = false;
  static final String CSV_NAME = "Rooftop_FE_Mode" + MODE + ".csv";

  public static void main(String[] args) { run(); }

  public static Model run() {
    model = ModelUtil.create("Model");
    String outDir = findOutDir();
    if (outDir != null) { try { model.modelPath(outDir); } catch (Exception e) { } say("Output folder: " + outDir); }
    else say("No output folder found; results go to Results > Tables > FE results.");
    say("MODE " + MODE + (QUICK ? " (QUICK test)" : ""));

    List<String> rows = new ArrayList<String>();
    rows.add("mode,alpha,lc_m,a_over_lc,delta_over_lc,dx_over_lc,deltaMid_over_lc,flat,"
           + "U1,U2,G_energy_over_Gth,Gm_closure_over_Gth,Gt_closure_over_Gth,"
           + "G_closure_over_Gth,closure_minus_energy_rel,domain");

    try {
      buildParameters();
      buildGeometry();
      buildSelections();
      buildPhysics();
      buildMesh();
      buildStudiesAndCouplings();

      double[] aList = QUICK ? new double[]{10} : A_LIST;
      for (double aL : aList) {
        A_OVER_LC = aL;
        double[] ds = QUICK ? new double[]{1.0} : (aL == 10 ? DELTAS_A10 : DELTAS_A100);
        rows.add(runCase(ALPHA, 1.0, true));
        writeCsv(outDir, rows);
        for (double d : ds) {
          rows.add(runCase(ALPHA, d, false));
          writeCsv(outDir, rows);
        }
      }
    } catch (Exception e) {
      say("STOPPED WITH ERROR: " + firstLine(e));
      String[] lines = String.valueOf(e.getMessage()).split("\n");
      for (int i = 1; i < Math.min(lines.length, 12); i++) say("   | " + lines[i].trim());
    }

    if (outDir != null) {
      try { model.save(outDir + File.separator + "Rooftop_FE_Mode" + MODE + ".mph"); say("Saved .mph (last case)."); }
      catch (Exception e) { say("Could not save .mph: " + firstLine(e)); }
      try {
        PrintWriter w = new PrintWriter(new FileWriter(outDir + File.separator + "Rooftop_FE_Mode" + MODE + "_log.txt"));
        for (String s : logBuf) w.println(s);
        w.close();
      } catch (Exception e) { }
    }
    flushLog();
    return model;
  }

  // =====================================================================
  static void buildParameters() {
    model.param().set("E0", "200[GPa]");
    model.param().set("nu0", "0.3");
    model.param().set("h0", "1[mm]");
    model.param().set("alpha0", "0.02", "Facet slope alpha (sets lc)");
    model.param().set("alpha_g", "alpha0", "Slope used in the geometry (0 = flat calibration)");
    model.param().set("lc", "h0/(alpha0*sqrt(3*(1-nu0^2)))", "Fold length lc");
    model.param().set("a_over_lc", String.valueOf(A_OVER_LC));
    model.param().set("a_crack", "a_over_lc*lc", "Crack half-length a");
    model.param().set("delta_over_lc", "1");
    model.param().set("delta_tip", "delta_over_lc*lc", "Near tip to ridge distance");
    model.param().set("dx_seg", "min(lc/20, delta_tip/5)", "Released segment length");
    model.param().set("W0", DOMAIN + "*a_crack");
    model.param().set("L0", DOMAIN + "*a_crack");
    model.param().set("ridge_len", "min(5*a_crack, 50*lc)", "Refined length of the ridge");
    model.param().set("S_inf", SINF + "[N/m]", "Remote in-plane shear (mode II)");
    model.param().set("V0", V0 + "[N/m]", "Face transverse shear (mode III)");
    model.param().set("D0", "E0*h0^3/(12*(1-nu0^2))");
  }

  // =====================================================================
  static void buildGeometry() {
    model.component().create("comp1", true);
    model.component("comp1").geom().create("geom1", 3);
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
    model.component("comp1").geom("geom1").create("ext1", "Extrude");
    model.component("comp1").geom("geom1").feature("ext1").set("workplane", "wp1");
    model.component("comp1").geom("geom1").feature("ext1").selection("input").set(new String[]{"wp1"});
    model.component("comp1").geom("geom1").feature("ext1").setIndex("distance", "ridge_len", 0);
    model.component("comp1").geom("geom1").feature("ext1").setIndex("distance", "L0", 1);
    model.component("comp1").geom("geom1").feature("ext1").set("reverse", true);
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
  static final String[] BOX1 = {"selFar", "selCrack", "selSeg", "selNear", "selRight", "selTop", "selRidge", "selSideL", "selSideR"};
  static final String[] BOX0 = {"selTipPts", "selOrigin", "selFarPt", "selRidgeTop"};
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
    double t = 0.01 * dx, tL = 1e-6 * L0;
    setBox("selFar",   -W0 - t, c1 + t,       -t, t);
    setBox("selCrack",  c1 - t, c2 + t,       -t, t);
    setBox("selSeg",    c2 - t, c2 + dx + t,  -t, t);
    setBox("selNear",   c2 + dx - t, t,       -t, t);
    setBox("selRight", -t, W0 + t,            -t, t);
    setBox("selTop",   -W0 - tL, W0 + tL,     L0 - tL, L0 + tL);
    setBox("selSideL", -W0 - tL, -W0 + tL,    -tL, L0 + tL);
    setBox("selSideR",  W0 - tL,  W0 + tL,    -tL, L0 + tL);
    setBox("selRidge", -t, t,                 -t, Math.min(5 * a, 50 * lc) + tL);
    setBox("selTipPts", c2 - t, c2 + dx + t,  -t, t);
    setBox("selOrigin", -t, t,                -t, t);
    setBox("selFarPt",  W0 - tL, W0 + tL,     -t, t);
    setBox("selRidgeTop", -tL, tL,            L0 - tL, L0 + tL);
    setBox("selFaces", -W0 - tL, W0 + tL,     -tL, L0 + tL);
  }

  // =====================================================================
  // Facet direction t = (1, 0, s)/q and normal n = (-s, 0, 1)/q, with
  // s = alpha_g*sign(-x) (s = +alpha on the left facet, -alpha on the right)
  // and q = sqrt(1+s^2).
  static final String S_EXPR = "(alpha_g*sign(-x))";
  static final String Q_EXPR = "sqrt(1+alpha_g^2)";

  static void buildPhysics() {
    model.component("comp1").material().create("mat1", "Common");
    model.component("comp1").material("mat1").propertyGroup("def").set("youngsmodulus", "E0");
    model.component("comp1").material("mat1").propertyGroup("def").set("poissonsratio", "nu0");
    model.component("comp1").material("mat1").propertyGroup("def").set("density", "1000[kg/m^3]");

    model.component("comp1").physics().create("shell", "Shell", "geom1");
    model.component("comp1").physics("shell").feature("to1").set("d", "h0");

    if (MODE.equals("II")) {
      // top edge: S_inf along t
      model.component("comp1").physics("shell").create("elTop", "EdgeLoad", 1);
      model.component("comp1").physics("shell").feature("elTop").selection().named("selTop");
      model.component("comp1").physics("shell").feature("elTop").set("forceReferenceLength",
        new String[]{"S_inf/" + Q_EXPR, "0", "S_inf*" + S_EXPR + "/" + Q_EXPR});
      // side edges: +S_inf along y on x = +W0, -S_inf on x = -W0
      model.component("comp1").physics("shell").create("elR", "EdgeLoad", 1);
      model.component("comp1").physics("shell").feature("elR").selection().named("selSideR");
      model.component("comp1").physics("shell").feature("elR").set("forceReferenceLength", new String[]{"0", "S_inf", "0"});
      model.component("comp1").physics("shell").create("elL", "EdgeLoad", 1);
      model.component("comp1").physics("shell").feature("elL").selection().named("selSideL");
      model.component("comp1").physics("shell").feature("elL").set("forceReferenceLength", new String[]{"0", "-S_inf", "0"});
    } else {
      // crack face y = 0+: V0 along the facet normal n (left facet, s = +alpha_g)
      model.component("comp1").physics("shell").create("elFace", "EdgeLoad", 1);
      model.component("comp1").physics("shell").feature("elFace").selection().named("selCrack");
      model.component("comp1").physics("shell").feature("elFace").set("forceReferenceLength",
        new String[]{"-V0*alpha_g/" + Q_EXPR, "0", "V0/" + Q_EXPR});
    }

    addAnti("pcU", "pcW", "pcR", "selSym");
    addAnti("pcSegU", "pcSegW", "pcSegR", "selSeg");

    // rigid-body supports (carry no force)
    model.component("comp1").physics("shell").create("disp1", "Displacement0", 0);
    model.component("comp1").physics("shell").feature("disp1").selection().named("selOrigin");
    model.component("comp1").physics("shell").feature("disp1").set("Direction", new String[]{"free", "prescribed", "free"});
    model.component("comp1").physics("shell").create("disp2", "Displacement0", 0);
    model.component("comp1").physics("shell").feature("disp2").selection().named("selFarPt");
    model.component("comp1").physics("shell").feature("disp2").set("Direction", new String[]{"free", "prescribed", "free"});
    model.component("comp1").physics("shell").create("disp3", "Displacement0", 0);
    model.component("comp1").physics("shell").feature("disp3").selection().named("selRidgeTop");
    model.component("comp1").physics("shell").feature("disp3").set("Direction", new String[]{"free", "free", "prescribed"});
  }

  static void addAnti(String tU, String tW, String tR, String sel) {
    String[][] c = {{tU, "u"}, {tW, "w"}, {tR, "arx+" + S_EXPR + "*arz"}};
    for (String[] k : c) {
      model.component("comp1").physics("shell").create(k[0], "PointwiseConstraint", 1);
      model.component("comp1").physics("shell").feature(k[0]).selection().named(sel);
      model.component("comp1").physics("shell").feature(k[0]).set("constraintExpression", k[1]);
    }
    say("Antisymmetry via pointwise constraints u = 0, w = 0, ar.t = 0 on " + sel);
  }

  static void setSegActive(boolean on) {
    for (String t : new String[]{"pcSegU", "pcSegW", "pcSegR"})
      model.component("comp1").physics("shell").feature(t).active(on);
  }

  // =====================================================================
  static void buildMesh() {
    model.component("comp1").mesh().create("mesh1");
    model.component("comp1").mesh("mesh1").feature("size").set("custom", "on");
    String[][] local = {{"szTip", "0", "selTipPts"}, {"szSeg", "1", "selSeg"}, {"szNear", "1", "selNear"},
                        {"szCrack", "1", "selCrack"}, {"szRidge", "1", "selRidge"}};
    for (String[] s : local) {
      model.component("comp1").mesh("mesh1").create(s[0], "Size");
      model.component("comp1").mesh("mesh1").feature(s[0]).selection().geom("geom1", Integer.parseInt(s[1]));
      model.component("comp1").mesh("mesh1").feature(s[0]).selection().named(s[2]);
      model.component("comp1").mesh("mesh1").feature(s[0]).set("custom", "on");
      model.component("comp1").mesh("mesh1").feature(s[0]).set("hmaxactive", true);
      model.component("comp1").mesh("mesh1").feature(s[0]).set("hgradactive", true);
      model.component("comp1").mesh("mesh1").feature(s[0]).set("hgrad", "1.15");
    }
    model.component("comp1").mesh("mesh1").create("ftri1", "FreeTri");
    model.component("comp1").mesh("mesh1").feature("ftri1").selection().geom("geom1", 2);
    model.component("comp1").mesh("mesh1").feature("ftri1").selection().named("selFaces");
  }

  static void updateMesh(double lc, double a, double delta, double dx) {
    model.component("comp1").mesh("mesh1").feature("size").set("hmax", String.valueOf(a));
    model.component("comp1").mesh("mesh1").feature("size").set("hmin", String.valueOf(dx / 20));
    model.component("comp1").mesh("mesh1").feature("size").set("hgrad", "1.2");
    model.component("comp1").mesh("mesh1").feature("szTip").set("hmax",   String.valueOf(dx / 4));
    model.component("comp1").mesh("mesh1").feature("szSeg").set("hmax",   String.valueOf(dx / 4));
    model.component("comp1").mesh("mesh1").feature("szNear").set("hmax",  String.valueOf(Math.min(lc / 4, Math.max(dx / 4, (delta - dx) / 4))));
    model.component("comp1").mesh("mesh1").feature("szCrack").set("hmax", String.valueOf(a / 20));
    model.component("comp1").mesh("mesh1").feature("szRidge").set("hmax", String.valueOf(lc / 2));
    model.component("comp1").mesh("mesh1").run();
  }

  // =====================================================================
  static void buildStudiesAndCouplings() {
    String[][] ops = {{"intTop", "selTop"}, {"intSL", "selSideL"}, {"intSR", "selSideR"}, {"intCrack", "selCrack"}};
    for (String[] o : ops) {
      model.component("comp1").cpl().create(o[0], "Integration");
      model.component("comp1").cpl(o[0]).selection().geom("geom1", 1);
      model.component("comp1").cpl(o[0]).selection().named(o[1]);
      model.component("comp1").cpl(o[0]).set("opname", o[0]);
    }
    model.component("comp1").cpl().create("sumSeg", "Integration");
    model.component("comp1").cpl("sumSeg").selection().geom("geom1", 1);
    model.component("comp1").cpl("sumSeg").selection().named("selSeg");
    model.component("comp1").cpl("sumSeg").set("opname", "sumSeg");
    model.component("comp1").cpl("sumSeg").set("method", "summation");
    model.study().create("std1"); model.study("std1").create("stat", "Stationary");
    model.study().create("std2"); model.study("std2").create("stat", "Stationary");
  }

  // U of the half model = 1/2 * work of the applied loads
  static double energy(String dset) {
    if (MODE.equals("II"))
      return eval(dset, "0.5*(intTop(S_inf*(u+" + S_EXPR + "*w)/" + Q_EXPR + ") + intSR(S_inf*v) - intSL(S_inf*v))");
    return eval(dset, "0.5*intCrack(V0*(-alpha_g*u+w)/" + Q_EXPR + ")");
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
    try {
      model.sol(sol).feature("s1").feature("dDef").set("linsolver", "pardiso");
      try { model.sol(sol).feature("s1").feature("dDef").set("errorchk", "off"); } catch (Exception e2) { }
      model.sol(sol).feature("s1").feature("fc1").set("linsolver", "dDef");
    } catch (Exception e) { say("   (could not force the direct solver: " + firstLine(e) + ")"); }
    model.sol(sol).runAll();
    try { model.result().dataset().remove(dset); } catch (Exception e) { }
    model.result().dataset().create(dset, "Solution");
    model.result().dataset(dset).set("solution", sol);
  }

  static double eval(String dset, String expr) {
    String t = "gevTmp";
    try { model.result().numerical().remove(t); } catch (Exception e) { }
    model.result().numerical().create(t, "EvalGlobal");
    model.result().numerical(t).set("data", dset);
    model.result().numerical(t).set("expr", new String[]{expr});
    return model.result().numerical(t).getReal()[0][0];
  }

  // =====================================================================
  static String runCase(double alpha, double dOverLc, boolean flat) {
    double lc = H / (alpha * Math.sqrt(3 * (1 - NU * NU)));
    double a = A_OVER_LC * lc, delta = dOverLc * lc;
    double dx = Math.min(lc / 20, delta / 5);
    double D = E * H * H * H / (12 * (1 - NU * NU)), c3 = (1 - NU) * (3 + NU);
    double Gth = MODE.equals("II") ? Math.PI * SINF * SINF * a / (E * H * H)
                                   : Math.PI * V0 * V0 * a * a * a / (4 * c3 * D * H);

    say(String.format("---- mode %s  a/lc=%.4g  delta/lc=%.4g  %s  (lc=%.4g m, dx/lc=%.3g)",
        MODE, A_OVER_LC, dOverLc, flat ? "FLAT CALIBRATION" : "", lc, dx / lc));
    model.param().set("a_over_lc", String.valueOf(A_OVER_LC));
    model.param().set("alpha0", String.valueOf(alpha));
    model.param().set("alpha_g", flat ? "0" : "alpha0");
    model.param().set("delta_over_lc", String.valueOf(dOverLc));

    model.component("comp1").geom("geom1").run();
    updateSelections(lc, a, delta, dx);
    checkSelections();
    updateMesh(lc, a, delta, dx);
    meshReport();

    setSegActive(false);
    solve("std2", "sol2", "dset2");
    setSegActive(true);
    solve("std1", "sol1", "dset1");
    solsExist = true;

    double U1 = energy("dset1"), U2 = energy("dset2");
    double Ge = 2 * (U2 - U1) / (H * dx);

    // Closure split. On the segment x < 0, so s = alpha_g.
    String q = Q_EXPR;
    String Rt = "(reacf(u)+alpha_g*reacf(w))/" + q, dt = "(withsol('sol2',u)+alpha_g*withsol('sol2',w))/" + q;
    String Rn = "(-alpha_g*reacf(u)+reacf(w))/" + q, dn = "(-alpha_g*withsol('sol2',u)+withsol('sol2',w))/" + q;
    String rot = "reacf(arx)*withsol('sol2',arx)+reacf(ary)*withsol('sol2',ary)+reacf(arz)*withsol('sol2',arz)";
    double Gm = Double.NaN, Gt = Double.NaN;
    try {
      Gm = -eval("dset1", "sumSeg(" + Rt + "*" + dt + ")") / (H * dx);
      Gt = -eval("dset1", "sumSeg(" + Rn + "*" + dn + "+" + rot + ")") / (H * dx);
    } catch (Exception e) { say("Closure evaluation failed: " + firstLine(e)); }
    // As in mode I: in the half model the reaction is the full force between
    // the faces and the displacement is half the jump, so no extra factor.
    // Closure must equal energy (both use the same two states).
    double Gc = Gm + Gt, rel = (Gc - Ge) / Ge;
    if (U1 == U2) say("WARNING: U1 == U2. The segment release did not change the model.");
    say(String.format("   G_energy/G_th = %.5f   Gm/G_th = %.5f   Gt/G_th = %.5f   (closure-energy)/energy = %.2e",
        Ge / Gth, Gm / Gth, Gt / Gth, rel));

    double dMid = (delta - dx / 2) / lc;
    return String.format(Locale.US, "%s,%.4g,%.6e,%.4g,%.6g,%.6g,%.6g,%d,%.10e,%.10e,%.6f,%.6f,%.6f,%.6f,%.3e,%.4g",
        MODE, alpha, lc, A_OVER_LC, dOverLc, dx / lc, dMid, flat ? 1 : 0,
        U1, U2, Ge / Gth, Gm / Gth, Gt / Gth, Gc / Gth, rel, DOMAIN);
  }

  static void checkSelections() {
    String[] all = {"selFar", "selCrack", "selSeg", "selNear", "selRight", "selTop", "selRidge", "selSideL", "selSideR",
                    "selTipPts", "selOrigin", "selFarPt", "selRidgeTop", "selFaces"};
    int[] expected = {1, 1, 1, 1, 1, 3, 1, 2, 2, 2, 1, 1, 1, 6};
    int[] dims     = {1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 2};
    for (int i = 0; i < all.length; i++) {
      int n;
      try { n = model.component("comp1").selection(all[i]).entities(dims[i]).length; }
      catch (Exception e) { say("(selection check skipped: " + e.getMessage() + ")"); return; }
      if (n != expected[i]) say("WARNING: selection " + all[i] + " has " + n + " entities, expected " + expected[i]);
      if (n == 0) throw new RuntimeException("Empty selection " + all[i]);
    }
  }

  static void writeCsv(String dir, List<String> rows) {
    updateTable(rows);
    if (rows.size() == 2) say("CSV " + rows.get(0));
    say("CSV " + rows.get(rows.size() - 1));
    if (dir == null) return;
    try {
      PrintWriter w = new PrintWriter(new FileWriter(dir + File.separator + CSV_NAME));
      for (String r : rows) w.println(r);
      w.close();
    } catch (Exception e) { say("CSV write failed: " + e.getMessage()); }
  }

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

  static String firstLine(Exception e) {
    String m = String.valueOf(e.getMessage());
    int i = m.indexOf('\n');
    return i > 0 ? m.substring(0, i) : m;
  }

  static void meshReport() {
    int nt = -1;
    try { nt = model.component("comp1").mesh("mesh1").getNumElem("tri"); say("   mesh: " + nt + " triangles"); }
    catch (Exception e) { say("   (mesh report unavailable: " + firstLine(e) + ")"); }
    if (nt == 0) throw new RuntimeException("The mesh is empty.");
  }

  static boolean tableMade = false;
  static void updateTable(List<String> rows) {
    try {
      if (!tableMade) { model.result().table().create("tblFE", "Table"); model.result().table("tblFE").label("FE results"); tableMade = true; }
      int nr = rows.size() - 1;
      String[] head = rows.get(0).split(",");
      double[][] data = new double[nr][head.length];
      for (int i = 0; i < nr; i++) {
        String[] f = rows.get(i + 1).split(",");
        for (int j = 0; j < head.length; j++) { try { data[i][j] = Double.parseDouble(f[j]); } catch (Exception e) { data[i][j] = Double.NaN; } }
      }
      model.result().table("tblFE").setColumnHeaders(head);
      model.result().table("tblFE").setTableData(data);
    } catch (Exception e) { say("Table update failed: " + e.getMessage()); }
  }

  static final String[] OUT_DIRS = {
    "/path/to/folded-sheet-fracture/results/comsol",
    "/path/to/folded-sheet-fracture/results/comsol",
    "../../results/comsol",        // relative to scripts/comsol, where the batch is started
    "."
  };
  static String findOutDir() {
    for (String d : OUT_DIRS) { try { if (new File(d).isDirectory()) return d; } catch (Exception e) { } }
    return null;
  }
}

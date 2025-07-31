(* ::Package:: *)

(* TunnelingActions Package *)

(* Implementation of arxiv 1805.03680 *)

BeginPackage["TunnelingActions`"]

Vt3::usage           = "Vt3[Vfunc_, phi_, phi0_, d_] gives the cubic approximation of the tunneling potential.";
Vt4::usage           = "Vt4[a4_, phi_, phi0_] gives the quartic correction term.";
Vt::usage            = "Vt[Vfunc_, a4_, phi_, phi0_, d_] combines Vt3 and Vt4 for the full approximate potential.";
Coefficients::usage  = "Coefficients[Vfunc_, phiT_, phi0_, d_] returns the {ad, bd, cd} coefficients of the quartic equation.";
SolveA4::usage       = "SolveA4[Vfunc_, phiT_, phi0_, d_] solves for the quartic coupling a4.";
Action::usage        = "Action[Vfunc_, phi0_, d_] approximates S_d via a 1D phi-integral.";
BubbleProfile::usage = "BubbleProfile[Vfunc_, phi0_, d_, dr_:0.001] returns {(r,phi)} via the thin-wall formula.";

Begin["`Private`"]

(* Cubic term around phi0 *)
Vt3[Vfunc_, phi_, phi0_, d_] := Module[{V0, Vp0},
  V0 = Vfunc[phi0]; Vp0 = D[Vfunc[x], x] /. x -> phi0;
  V0/phi0*phi +
  phi/(d phi0^2)*((d-1)*phi0 Vp0 - d V0)*(phi - phi0) +
  phi/(d phi0^3)*((d-1)*phi0 Vp0 - 2 d V0)*(phi - phi0)^2
];

(* Quartic correction term *)
Vt4[a4_, phi_, phi0_] := a4 phi^2 (phi - phi0)^2;

(* Full quartic approx. potential *)
Vt[Vfunc_, a4_, phi_, phi0_, d_] := Vt3[Vfunc, phi, phi0, d] + Vt4[a4, phi, phi0];

(* Auxiliary Q-derivatives at phiT *)
Q1[phiT_, phi0_] := 2 phiT (phiT - phi0)^2 + 2 phiT^2 (phiT - phi0);
Q2[phiT_, phi0_] := 2 (phiT - phi0)^2 + 8 phiT (phiT - phi0) + 2 phiT^2;

(* Compute coefficients ad, bd, cd *)
Coefficients[Vfunc_, phiT_, phi0_, d_] := Module[
  {VT, v3, vp, vpp, q, q1, q2},
  VT = Vfunc[phiT] - Vfunc[0];
  v3 = Vt3[Vfunc, phiT, phi0, d];
  vp = D[Vt3[Vfunc, x, phi0, d], x] /. x -> phiT;
  vpp = D[Vt3[Vfunc, x, phi0, d], {x,2}] /. x -> phiT;
  q = phiT^2 (phiT - phi0)^2;
  q1 = Q1[phiT, phi0]; q2 = Q2[phiT, phi0];
  {3 q1^2 - 4 q q2,
   -2 vpp q + 3 vp q1 - 2 q2 (v3 - VT),
   3 vp^2 - 4 vpp (v3 - VT)}
];

(* Solve for quartic coefficient a4 *)
SolveA4[Vfunc_, phiT_, phi0_, d_] := Module[{ad, bd, cd, disc},
  {ad, bd, cd} = Coefficients[Vfunc, phiT, phi0, d];
  disc = bd^2 - ad cd;
  If[disc < 0, Indeterminate, (-bd - Sqrt[disc])/ad]
];

(* 1D phi-integral action approximation *)
Action[Vfunc_, phi0_, d_] := Module[
  {phiT, a4, vt, vtprime, prefac},
  phiT = FindMaximum[Vfunc[x], {x, phi0/2}][[2, 1, 2]];
  a4 = SolveA4[Vfunc, phiT, phi0, d];
  vt[phi_] = Vt[Vfunc, a4, phi, phi0, d];
  vtprime[phi_] = D[vt[phi], phi];
  prefac = ((d - 1)^(d - 1) (2 Pi)^(d/2))/Gamma[d/2 + 1];
  prefac NIntegrate[
    (Vfunc[phi] - vt[phi])^(d/2)/Abs[vtprime[phi]]^(d - 1),
    {phi, 0, phi0}, Method -> "GlobalAdaptive"]
];

(* Thin-wall bounce profile (r vs phi) *)
BubbleProfile[Vfunc_, phi0_, d_, dr_:0.001] := Module[
  {phiT, a4, vt, vtprime, phiList, rList},
  phiT = FindMaximum[Vfunc[x], {x, phi0/2}][[2, 1, 2]];
  a4 = SolveA4[Vfunc, phiT, phi0, d];
  vt[phi_] = Vt[Vfunc, a4, phi, phi0, d];
  vtprime[phi_] = D[vt[phi], phi];
  phiList = Range[dr, phi0, dr];
  rList = 3 Table[Sqrt[2 (Vfunc[ph] - vt[ph])/vtprime[ph]^2], {ph, phiList}];
  Transpose[{rList, phiList}]
];

End[]
EndPackage[]

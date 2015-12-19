BeginPackage["StoppingPower`"]

dedxc::usage="dedxc[vp,(pg:Automatic)]"
dedtc::usage="dedtc[vp,(pg:Automatic)]"
dedxci::usage="dedxci[vp,(pg:Automatic)]"
dedtci::usage="dedxci[vp,(pg:Automatic)]"

dedxqm::usage="dedxqm[vp]"
dedtqm::usage="dedtqm[vp]"
dedxqmi::usage="dedxqmi[vp]"
dedtqmi::usage="dedxqmi[vp]"


dedx::usage="dedx[vp,(pg:Automatic),(umax:Infinity)] 
To change PrecisionGoal only (for CL): dedx[vp,pg] ; 
To change umax only (for QM): dedx[vp,Automatic,umax] ;
To change PrecisionGoal and umax: dedx[vp,pg,umax] "
dedt::usage="dedt[vp,(pg:Automatic),(umax:Infinity)] 
To change PrecisionGoal only (for CL): dedt[vp,pg] ; 
To change umax only (for QM): dedt[vp,Automatic,umax] ;
To change PrecisionGoal and umax: dedt[vp,pg,umax] "
dedxi::usage"dedxi[vp]"
dedti::usage"dedti[vp]"

param::usage="param[\[Beta], Zp, mp, mb, Zb, nb, K, r, i]"

dedxLP::usage="dedxLP[vp]"
dedtLP::usage="dedxLP[vp]"

dedxiLP::usage="dedxLPb[vp]"
dedtiLP::usage="dedxLPb[vp]"

paramLP::usage="paramLP[\[Beta], Zp, mp, mb, Zb, nb, K, r, i,g0,h0]"

(* Debugging *)


coulombLogLP::usage=""
daw::usage=""
dawb::usage=""
Fr::usage=""
Fi::usage=""
Farg::usage=""
Fabs::usage=""
H::usage=""
Hb::usage=""
inttwo::usage=""
inttwob::usage=""
intone::usage=""
dint1::usage=""
dint2::usage=""
dint3::usage=""
dint4::usage=""

repsi::usage=""
dintQM1::usage=""
dintQM2::usage=""
dintQM::usage=""
nintQM1::usage=""
nintQM2::usage=""
intQM1::usage=""
intQM2::usage=""
intQM::usage=""


Begin["`PrivateBPS`"]
(* Brown-Preston-Singleton: dE/dx *)

(* Clean the Slate *)
Clear[daw,dawb,Fr,Fi,Fabs,Farg,H, inttwo]
Clear[dintA,dintB,dintC,dintD,intone]
Clear[dedxc,dedtc,repsi,nintQM1,nintQM2]
Clear[intQM1,intQM2,intQM,dedxqm,dedtqm]
Clear[dedx,dedt,dedxci,dedtci,dedxqmi,dedtqmi]
Clear[dedxi, dedti, dedxiLP, dedtiLP]

(* Classical: *)

  sqrt = Sqrt[Pi] // N;
  sqrt2= sqrt/2   // N;

  (* note: daw[x]=NIntegrate[Exp[y^2 - x^2], {y, 0, x}]  *) 

  daw[x_]   :=  sqrt2*Exp[-x^2]*Erfi[x]  (* R -> R       *)
  dawb[xx_] := Map[daw, xx]              (* List -> List *)

  Fr[x_, kb2_] := Plus @@ (kb2*(1 - 2*x*dawb[x]))
  Fi[x_, kb2_] := sqrt*Plus @@ (kb2*x*Exp[-x^2])
  Fabs[x_, kb2_] := Sqrt[Fr[x, kb2]^2 + Fi[x, kb2]^2]
  Farg[x_, kb2_] := Arg[Fr[x, kb2] + I*Fi[x, kb2]]
  H[abv_, kb2_, K_] := 
    -2*(Fi[abv, kb2]*Log[Fabs[abv, kb2]/K^2] + 
    Fr[abv, kb2]*Farg[abv, kb2])

  inttwo[abv_, kb2_, K_, pg_:Automatic,ZERO_:0.000001]:= 
    2*NIntegrate[u H[abv u, kb2, K], {u, ZERO,1}, 
		 PrecisionGoal -> pg] (*somthing wrong here*)

  Hb[abv_, kb2_, K_] := sqrt*kb2*abv*Exp[-abv*abv]*H[abv,kb2,K]/Fi[abv,kb2]
  inttwob[abv_, kb2_, K_, pg_:Automatic,ZERO_:0.000001]:= 
    Table[
      2*NIntegrate[u Hb[abv u, kb2, K][[i]], {u, ZERO,1}, 
      PrecisionGoal -> pg],
     {i,1,Nb}
  ]

  dint1[a2_,u_] := Exp[-a2*u]*Log[1-u]/Sqrt[u]
  dint2[a2_,u_] := Exp[-a2*u]*Log[1-u]*Sqrt[u]
  dint4[a2_,u_] := Exp[-a2*u]Log[u]*Sqrt[u] 
  dint3[a2_,u_] := (Exp[-a2*u]-1+a2*u)Log[u]/Sqrt[u] 
  (* dint3 is non-singular at u=0, and the (-1+a2*u) piece integrates to -4+4*a2/9 *)
  intone[a_, b_, c_, pg_:Automatic] :=
  Module[{a2,iABCD,intex,int1,int2,int3,int4,umax,amax},
	a2=a^2;
        amax=50;
        int1=If[a<amax,Sqrt[Pi]*(-((EulerGamma*Erf[a])/a) +                     \
             Derivative[0, 1, 0][Hypergeometric1F1Regularized][1/2, 3/2, -a2]), \
             -Sqrt[Pi]/(2*a^3)];
        int2=If[a<amax, EulerGamma/(a2*E^a2) - (EulerGamma*Sqrt[Pi]*Erf[a])/(2*a^3) +          \
          (Sqrt[Pi]*Derivative[0, 1, 0][Hypergeometric1F1Regularized][1/2, 3/2, -a2])/(2*a2) - \
          (Sqrt[Pi]*Derivative[0, 1, 0][Hypergeometric1F1Regularized][3/2, 3/2, -a2])/(2*a2),  \
          -3*Sqrt[Pi]/(4*a^5)];
	int3=If[a<amax, NIntegrate[dint3[a2,u],{u,0,1},PrecisionGoal -> 10], \
             -Sqrt[Pi](Log[4*a2] + EulerGamma)/a - (-4+4*a2/9)];
        int4=If[a<amax, NIntegrate[dint4[a2,u],{u,0,1},PrecisionGoal -> pg], \
               -Sqrt[Pi](Log[4*a2] + EulerGamma - 2)/(2*a^3)];
        iABCD=int3-int1+b*(int2-int4);
	intex=(2-c+b*c/(2*a2))*Sqrt[Pi]*Erf[a]/a-b*c*Exp[-a2]/a2+(-4+4*a2/9); 
	iABCD+intex
    ]

  (* vp is measured in units of the thermal velocity vth *)
  dedxc[vp_, pg_:Automatic] := Module[{abv, bbv, tint1, i1}, 
      abv = ab*vp;
      bbv = bb*(vp^2);
      tint1 = Table[intone[abv[[j]],bbv[[j]],cb[[j]],pg],{j,1,Nb}];
      i1 = Plus @@ ((cp1/vp)*kb2*eb*tint1);
      i1 = i1 + cp2*inttwo[abv, kb2, K, pg];
      i1 = i1 - (cp3/vp^2)*H[abv,kb2,K];
      Re[i1*System`convfact]
  ]
  dedtc[vp_, pg_:Automatic] := dedxc[vp, pg]*vp*vth

  dedxci[vp_, pg_:Automatic] := Module[{abv, bbv, tint1, i1}, 
      abv = ab*vp;
      bbv = bb*(vp^2);
      i1 = Table[intone[abv[[j]],bbv[[j]],cb[[j]],pg],{j,1,Nb}];
      i1 = (cp1/vp)*kb2*eb*i1;
      i1 = i1 + cp2*inttwob[abv, kb2, K, pg];
      i1 = i1 - (cp3/vp^2)*Hb[abv,kb2,K];
      Re[i1*System`convfact]
  ];
   dedtci[vp_, pg_:Automatic] := dedxci[vp,pg]*vp*vth


(* Quantum:*)

  repsi[z_] := Re[PolyGamma[1 + I*z]]

  dintQM1[a_,e_,v_,u_] := Exp[-a*(v^2+u^2)]*(Log[Abs[e/u]] - repsi[e/u])*Sinh[2*a*v*u] 
  dintQM2[a_,e_,v_,u_] := Exp[-a*(v^2+u^2)]*(Log[Abs[e/u]] - 
   repsi[e/u])*(Cosh[2*a*v*u] - (Sinh[2*a*v*u]/(2*a*v*u)))/u


   dintQM[a_, e_, v_, u_]:= -rmb0*dintQM1[a,e,v,u] +
    rMb0*dintQM2[a,e,v,u]


  nintQM1[a_, e_, v_] := 
    Module[{eps,l,m},
	   eps=15./Sqrt[a];
	   m=v+eps;
           l=v-eps;
           l=Max[0,l];
        Re[(1/v)*NIntegrate[dintQM1[a,e,v,u], {u, l, m}]] 
      ]

  nintQM2[a_, e_, v_] := 
    Module[{eps,l,m},
	   eps=15./Sqrt[a];
	   m=v+eps;
           l=v-eps;
           l=Max[0,l];
           Re[NIntegrate[dintQM2[a,e,v,u], {u, l, m}]]
  ]

  intQM1[a_, e_, v_] :=
      Table[-rmb0[[j]]*nintQM1[a[[j]],e[[j]],v], 
      {j, 1, Length[a]}] //N

  intQM2[a_, e_, v_] :=
      Table[rMb0[[j]]*nintQM2[a[[j]],e[[j]],v], 
      {j, 1, Length[a]}] //N

  intQM[a_, e_, v_] := intQM1[a, e, v] + 
    intQM2[a, e, v]

  dedxqm[vp_] := (cp1*System`convfact/vp)*Plus @@ (kb2*fb*intQM[ab2, etb, vp])
  dedtqm[vp_] := dedxqm[vp]*vp*vth

  dedxqmi[vp_] := (cp1*System`convfact/vp)*(kb2*fb*intQM[ab2, etb, vp])
  dedtqmi[vp_] := dedxqmi[vp]*vp*vth 

(* Full = Classical + Quantum *)

  dedx[vp_, pg_:Automatic] := dedxc[vp, pg] + dedxqm[vp]
  dedt[vp_, pg_:Automatic] := dedx[vp, pg, ul, um]*vp*vth

  dedxi[vp_, pg_:Automatic] := dedxci[vp,pg] + dedxqmi[vp]
  dedti[vp_, pg_:Automatic] := dedxi[vp, pg]*vp*vth

(* Parameters and Factors *)

  System`KtoeV = 1/(38.68*300)    // N;
  System`cmtoa0 = 1/(5.3*10^(-9)) // N;
  System`mpgm = 1.67*10^(-27)*1000;     (* mass of proton   in gm *)
  System`mpev = 0.938271998*10^9;       (* mass of proton   in eV *)
  System`megm = 9.11*10^(-31)*1000;     (* mass of electron in gm *)
  System`meev = 0.510998902*10^6;       (* mass of electron in gm *)

  cmtoa0 = 1/(5.3*10^(-9)) // N;
  mtr = 10^(-6);                        (* distance unit in meters *)
  ev = 10^6;                            (* energy unit in eV *)
  System`convfact=cmtoa0*(mtr*100)/ev;  (* dE/dx from atomic units to MeV/\[Mu]m *)

  c = 2.998*10^(10);       (* speed-of-light in cm/s *)
  Be = 13.6;               (* binding energy of 1st Hydrogen electron in eV *)

  param[\[Beta]_, Zp_, mp_, mb_, Zb_, nb_, r_, i_, pm_] :=
    Module[{},

      Clear[kb2,kD,Nb,mm,mp0,mb0,Mpb0,mpb0,cp1,cp2,cp3,System`K];
      Clear[ab,bb,eb,cv,ab2,vthc,vth,etb,fb,rMb0,rmb0,gpb,gD];

      kb2 = 8*Pi*Be*\[Beta]*Zb*Zb*nb // N;  (* inverse Debey length squared:*)
      kD  = Sqrt[Plus @@ Abs[kb2]]   // N;  (* units a0^-2                  *)
      System`K = kD;           (* arbitrary wave number in units a_0^-1 *)

      Nb  = Length[mb];         (* number of plasma species        *)
      mm  = Part[mb,i] // N;    (* mass of index splama species    *)
      mp0 = mp/mm     // N;     (* rescale masses                  *)
      mb0 = mb/mm     // N;     (* rescale masses                  *)
      Mpb0= mp0 + mb0;          (* total mass ratio list           *)
      mpb0= mp0*mb0/Mpb0;       (* reduced mass ratio list         *)

      cp1 = (2*Be*Zp^2)          // N; (* parameter: units of eV - a0 *)
      cp2 = (Be*Zp^2)/(2*Pi)     // N; (* parameter: units of eV - a0 *)
      cp3 = (Be*Zp^2)/(r*Pi*mp0) // N; (* dimensionless parameter     *)

      ab = Sqrt[r*mb0/2]  // N;         (* dimensionless parameter lists *)
      bb = r*Mpb0         // N;         
      eb = (mb0/mp0)/Sqrt[2*Pi*r*mb0] // N;  
      cb = 2 - 2*EulerGamma - Log[Abs[(2*Be)*\[Beta]*Zp*Zb*K*mb0/mpb0]];
      ab2= ab*ab;

      vthc = Sqrt[r/(\[Beta]*mm)];         (* thermal velocity of mm: units of c *)
      vth  = c*vthc;                            (* thermal velocity of mm: units cm/s *)
      etb  = 2*Be*Abs[Zp*Zb]*(2.686/10^4)/vthc; (* quantum parameter *)
      fb   = 2/Sqrt[2*Pi*r*mb0] // N;
      rMb0 = Mpb0/mp0           // N; (* rm *)
      rmb0 = mb0/mp0            // N; (* rr *)
      If[pm==1, 
      gpb=Abs[2*Be*(\[Beta])*Zp*Zb*Sqrt[kb2]]//N;
      gD=Sqrt[Plus@@(gpb*gpb)]//N;
      Print["gD=",gD];
      Print["gpb=",gpb];    
      Print["etb*Sqrt[mb0]=",etb*Sqrt[mb0]];
     ];
  ]

End[]

Begin["`PrivateLP`"]
(* Li-Petrasso: dE/dx *)

  Clear[dedxiLP];
  Clear[dedtiLP];


  tspi = 2/Sqrt[Pi] // N;
  mu[x_] := Erf[Sqrt[x]] - tspi*Sqrt[x]*Exp[-x]
  dmu[x_] := tspi*Sqrt[x]*Exp[-x]
  mub[xx_List]  := Map[mu, xx]
  dmub[xx_List] := Map[dmu, xx]

  lnb[xx_List] := 
    Module[{logb, xxone},
      xxone = 1 + xx;
      logb = 1/(xxone*xxone) + qmb0/xxone;
      logb = 0.5*Log[logb] + Log[lb0];
      -logb
    ]
  coulombLogLP[vp_] := lnb[xb0*vp^2]

  gb[xx_List] := (1/System`r)*(lnb[xx]*(mub[xx]/mb0 - 
    dmub[xx]/mp0) + Erf[Sqrt[xx]]/mp0)

  temg = 2*Exp[-EulerGamma] // N;
  hb[xx_List] := UnitStep[(xx - 1)]*Log[temg*Sqrt[xx]]/(System`r*mb0)

  ghb[xx_List] := g0*gb[xx] + h0*hb[xx]

  dedxLP[vp_] := (cp1*convfact/vp^2)*(Plus @@ (kb2*ghb[xb0*vp^2]))
  dedtLP[vp_] := dedxLP[vp]*vp*vth

  dedxiLP[vp_] := (cp1*convfact/vp^2)*(kb2*ghb[xb0*vp^2])
  dedtiLP[vp_] := dedxiLP[vp]*vp*vth




(* Parameters and Factors *)


  System`KtoeV = 1/(38.68*300)    // N;
  System`cmtoa0 = 1/(5.3*10^(-9)) // N;
  System`mpgm = 1.67*10^(-27)*1000;     (* mass of proton   in gm *)
  System`mpev = 0.938271998*10^9;       (* mass of proton   in eV *)
  System`megm = 9.11*10^(-31)*1000;     (* mass of electron in gm *)
  System`meev = 0.510998902*10^6;       (* mass of electron in gm *)

  cmtoa0 = 1/(5.3*10^(-9)) // N;
  mtr = 10^(-6);                        (* distance unit in meters *)
  ev = 10^6;                            (* energy unit in eV *)
  System`convfact=cmtoa0*(mtr*100)/ev;  (* dE/dx from atomic units to MeV/\[Mu]m *)

  c = 2.998*10^(10);   (* speed-of-light in cm/s *)
  Be = 13.6;           (* binding energy of 1st Hydrogen electron in eV *)

  paramLP[\[Beta]_, Zp_, mp_, mb_List, Zb_, nb_, r_, i_,g_,h_] :=
    Module[{},
      cp1 = (2*Be*Zp^2)          // N; (* parameter: units of eV - a0 *)

      kb2 = 8*Pi*Be*\[Beta]*Zb*Zb*nb // N;  (* inverse Debye length squared: units a0^-2 *)
      kD  = Sqrt[Plus @@ Abs[kb2]]   // N;  (* total inverse Debye length: units a0^-1   *)

      Nb = Length[mb];          (* number of plasma species        *)
      mm = Part[mb,i] // N;     (* mass of index splama species    *)
      mp0 = mp/mm     // N;     (* rescale masses                  *)
      mb0 = mb/mm     // N;     (* rescale masses                  *)
      xb0 = r*mb0/2   // N;     (* rescaled x-variable             *)
      Mpb0 = mp0 + mb0;         (* total mass ratio list           *)
      mpb0 = mp0*mb0/Mpb0;      (* reduced mass ratio list         *)

      vthc = Sqrt[r/(\[Beta]*mm)];      (* thermal velocity of mm: units of c *)
      vth = c*vthc;                     (* thermal velocity of mm: units cm/s *)
      etb = 2*Be*Abs[Zp*Zb]*(2.686/10^4)/vthc; (* quantum parameter *)

      qmb0 = 1/(2*r*etb*etb*mb0) // N; 
      lb0 = Be*kD*Abs[Zp*Zb]*\[Beta]*(mb0/mpb0) // N; 
      g0=g;
      h0=h;
  ]

End[]

EndPackage[]

TCutG *myPionCentrCut_1998_uv(){
//========= Macro generated from object: PionCentrCut/
//========= by ROOT version6.28/06
   
   TCutG *cutg = new TCutG("CUTCentrPion_1998_uv",20);
   cutg->SetVarX("Eclusc_track");
   cutg->SetVarY("Dedx_track");
   cutg->SetTitle("");
   cutg->SetFillStyle(1000);
   cutg->SetPoint(0,0.0133238,7.96007);
   cutg->SetPoint(1,0.0369627,5.92014);
   cutg->SetPoint(2,0.0541547,4.14062);
   cutg->SetPoint(3,0.133668,3.27257);
   cutg->SetPoint(4,0.168052,2.83854);
   cutg->SetPoint(5,0.202436,2.75174);
   cutg->SetPoint(6,0.254011,2.66493);
   cutg->SetPoint(7,0.301289,2.53472);
   cutg->SetPoint(8,0.350716,2.31771);
   cutg->SetPoint(9,0.34212,2.31771);
   cutg->SetPoint(10,0.34212,2.10069);
   cutg->SetPoint(11,0.337822,1.97049);
   cutg->SetPoint(12,0.32063,1.44965);
   cutg->SetPoint(13,0.213181,1.36285);
   cutg->SetPoint(14,0.116476,1.57986);
   cutg->SetPoint(15,0.0412607,2.70833);
   cutg->SetPoint(16,0.00687678,1.05903);
   cutg->SetPoint(17,0.00902578,8.09028);
   cutg->SetPoint(18,0.0154728,8.09028);
   cutg->SetPoint(19,0.0133238,7.96007);
   
   return cutg;
}

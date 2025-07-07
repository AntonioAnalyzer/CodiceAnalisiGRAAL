TCutG *myProtonCentrCut_1999_uv(){
//========= Macro generated from object: myProtonCentrCut/Graph
//========= by ROOT version6.28/06
   
   TCutG *cutg = new TCutG("CUTCentrPro_1999_uv",25);
   cutg->SetVarX("Eclusc_track");
   cutg->SetVarY("Dedx_track");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->SetPoint(0,0.21533,3.48958);
   cutg->SetPoint(1,0.163754,4.31424);
   cutg->SetPoint(2,0.127221,5.05208);
   cutg->SetPoint(3,0.0863897,6.09375);
   cutg->SetPoint(4,0.0541547,7.82986);
   cutg->SetPoint(5,0.0326647,9.9566);
   cutg->SetPoint(6,0.0262178,11.3455);
   cutg->SetPoint(7,0.0219198,12.691);
   cutg->SetPoint(8,0.0197708,13.6892);
   cutg->SetPoint(9,0.0176218,15.0347);
   cutg->SetPoint(10,0.0111748,16.5538);
   cutg->SetPoint(11,0.0176218,17.4653);
   cutg->SetPoint(12,0.0391117,14.7743);
   cutg->SetPoint(13,0.0541547,11.9531);
   cutg->SetPoint(14,0.0885387,9.39236);
   cutg->SetPoint(15,0.125072,7.26562);
   cutg->SetPoint(16,0.180946,5.96354);
   cutg->SetPoint(17,0.266905,4.92187);
   cutg->SetPoint(18,0.34212,4.40104);
   cutg->SetPoint(19,0.40659,3.35937);
   cutg->SetPoint(20,0.333524,2.96875);
   cutg->SetPoint(21,0.27765,3.14236);
   cutg->SetPoint(22,0.247564,3.31597);
   cutg->SetPoint(23,0.213181,3.53299);
   cutg->SetPoint(24,0.21533,3.48958);
   
   return cutg;
}

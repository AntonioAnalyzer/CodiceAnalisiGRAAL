///////////////////////////////////////////////////
// Header file contenente l'argoritmo di processing
// l'algoritmo prende un albero e lo cerca di
// analizzare
///////////////////////////////////////////////////


void Analysis::process(TTree *alb, string root_file, string cartella){

  
  //Funzione di inizializzazione dell'albero passato
  Analysis::Init(alb);

  //Inizio algoritmo
  double RMP    = 0.9382720813;
  double RMN    = 0.9395654133;
  double RMD = 1.877;
  double CLIGHT = 29979245800.; //! cm/sec
  double DIST_WALL;
  double percent = 10.0;

  TCutG *ProtonCentrCut   = myProtonCentrCut();
  TCutG *PionCentrCut     = myPionCentrCut();
  TCutG *ProtonForwCut    = myProtonForwCut(cartella);
  TCutG *PionForwCut      = myPionForwCut();
  TCutG *DeuForwCut       = MyDeuForwCut(cartella);
  
  //TF1 *f = CutDeuForw(cartella);
  
  TH2D *protondxE         = new TH2D("protondxE","protondxE",100,0,1.2,100,0,20);
  TH2D *protondxEsoloPro  = new TH2D("protondxEsoloPro","protondxEsoloPro",100,0,1.2,100,0,20);
  TH2D *protondxEsoloPio  = new TH2D("protondxEsoloPio","protondxEsoloPio",100,0,1.2,100,0,20);
  TH2D *ddxtof_f          = new TH2D("ddxtof_f","ddxtof_f",100,0,0,100,0,0);
  TH2D *ddxtof_fow        = new TH2D("ddxtof_fow","ddxtof_fow",100,0,0,100,0,0);

  // questi servono nel caso si voglia sbloccare il codice sotto che produce una rootpla contenente i tagli in avanti delle particelle cariche dei diveris periodi

  string str1 = "protondxE_" + cartella;

  string str2 = "ddxtof_f_" + cartella;

  string str3 = "ddxtof_fow_" + cartella;

  string str4 = "Theta_trf_fch_" + cartella;
  string str5 = "Phi_trf_fch_" + cartella;
  string str6 = "Theta_trf_fn_" + cartella;
  string str7 = "Phi_trf_fn_" + cartella;
  string str8 = "N_show_" + cartella;
  string str9 = "E_show_"+ cartella;
  string str10 = "Iass_trf_" + cartella;
  string str11 = "ind_trf_" + cartella;
  string str12 = "../Forw_tracks_"+cartella;


  TH2D *protondxE         = new TH2D(str1.c_str(),str1.c_str(),100,0,1.2,100,0,20); // distribuzione perdita di energia vs energia nel BGO di tutte le tracce cariche(MWC-Bar-BGO) al centro
  TH2D *protondxEsoloPro  = new TH2D("protondxEsoloPro","protondxEsoloPro",100,0,1.2,100,0,20); // distribuzione perdita di energia vs energia nel BGO solo protoni dopo che questi sono stati identificati con il taglio sul grafico De/dx vs Er.
  TH2D *protondxEsoloPio  = new TH2D("protondxEsoloPio","protondxEsoloPio",100,0,1.2,100,0,20); // distribuzione perdita di energia vs energia nel BGO solo pioni dopo che questi sono stati identificati con il taglio sul grafico De/dx vs Er.
  TH2D *ddxtof_f          = new TH2D(str2.c_str(),str2.c_str(),100,0,0,100,0.,0.); // distribuzione perdita di energia vs tempo di volo per particelle neutre emesse in avanti.
  TH2D *ddxtof_fow        = new TH2D(str3.c_str(),str3.c_str(),50,0,0,50,0.,0.); // distribuzione perdita di energia vs tempo di volo per particelle cariche emesse in avanti.

  TH1D *N_show_ = new TH1D(str8.c_str(),"",100,0.,0.); // distribuzione del numero di colpi nel Russian Wall
  TH1D *E_show_ = new TH1D(str9.c_str(),"",20,0.,0.2); // distribuzione dell'energia depositata nel Russian Wall
  TH1D *Iass_trf_ = new TH1D(str10.c_str(),"",10,0.,10.); // distribuzione dell'indice che identifica il tipo di traccia in avanti
  TH1D *ind_trf = new TH1D(str11.c_str(),"",10,0.,10); // distribuzione dell'indice contenente l'informazione sull'associazione dei rivelatori in avanti
  TH1D *Theta_trf_fch = new TH1D(str4.c_str(),"",170,0.,0.); // distribuzione angolo theta delle tracce "cariche" emesse in avanti
  TH1D *Phi_trf_fch = new TH1D(str5.c_str(),"",350,0.,0.); // distribuzione angolo phi delle tracce "cariche" emesse in avanti
  TH1D *Theta_trf_fn = new TH1D(str6.c_str(),"",60,0.,0.); // distribuzione angolo theta delle tracce "cariche" emesse in avanti
  TH1D *Phi_trf_fn = new TH1D(str7.c_str(),"",72,0.,0.); // distribuzione angolo phi delle tracce "cariche" emesse in avanti

  TH3D *protondxE_Ind_bar = new TH3D("protondxE_Ind_bar","",32,0.5,32.5,200,0.,0.8,40,0.,40.); // grafico dE/dx vs Er vs Indice della barra	
  
  TFile *dati = new TFile(root_file.c_str(),"recreate");
  TTree *tree = new TTree("h80","Graal data Analysis");

  int idrun, idevt, ipol;

  vector<TLorentzVector> cphoton;
  vector<TLorentzVector> proton;
  // vector<TLorentzVector> cproton;
  vector<pair<double,double> > pionangles;
  vector<pair<double,double> > fphotonangles;
  // vector<TLorentzVector> fproton;
  vector<TLorentzVector> fneutron;
  vector<TLorentzVector> fDeuteron; // creiamo un vettore che contiene 
  //tutte le particelle promosse a "deuteroni in avanti" in ogni evento(ma lo facciamo solo se sulle analisi su deuterio)
  //vector<TLorentzVector> cDeuteron;  // lo stesso facciamo con i deuteroni al centro
  vector<TLorentzVector> ccharged; // ipotesi di un carico emesso ad angoli centrali che non sia un pione carico(da utilizzare in futuro..)
  vector<TLorentzVector> fcharged; // ipotesi di un carico emesso ad angoli in avanti che non sia un pione
  vector<pair<double,double> > fproangles;
  vector<pair<double,double> > fdeuangles;
  vector<double> Tof_fpro;
  vector<double> Tof_fdeu;
  vector<double> Tof_fch;
  int fch;

  
  TLorentzVector beam;

  tree->Branch("beam",    "TLorentzVector", &beam);
  tree->Branch("idrun",   &idrun, "idrun/I");   // numero del run
  tree->Branch("idevt",   &idevt, "idevt/I");   // numero dati dell'evento 
  tree->Branch("ipol",    &ipol,  "ipol/I");    // polarizzazione del fascio: 0 verticale, 1 orizzontale, 2 bremsstrahlung 
  tree->Branch("cphoton", &cphoton);
  // tree->Branch("cproton", &cproton);
  tree->Branch("proton",  &proton);
  // tree->Branch("fproton", &fproton);
  tree->Branch("fneutron",&fneutron);
  tree->Branch("pionangles",    &pionangles);
  tree->Branch("fphotonangles", &fphotonangles);
  tree->Branch("fDeuteron",&fDeuteron);

   // taglio energia de fascio vs angolo theta del protone per selezionare gli eventi del canale di fotoproduzione di etaprimo su protone sopra la soglia
   TCutG *cutg = new TCutG("cutg",22);
   cutg->SetVarX("Energy gamma(GeV)");
   cutg->SetVarY("Theta proton(degree)");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->SetPoint(0,1.44258,-0.171875);
   cutg->SetPoint(1,1.44196,2.69271);
   cutg->SetPoint(2,1.44917,5.46181);
   cutg->SetPoint(3,1.45606,7.2283);
   cutg->SetPoint(4,1.46515,8.99479);
   cutg->SetPoint(5,1.47988,11.0477);
   cutg->SetPoint(6,1.49586,13.1007);
   cutg->SetPoint(7,1.51529,15.3924);
   cutg->SetPoint(8,1.54256,17.5885);
   cutg->SetPoint(9,1.55979,19.0686);
   cutg->SetPoint(10,1.58048,20.3576);
   cutg->SetPoint(11,1.59928,21.3602);
   cutg->SetPoint(12,1.60336,21.0738);
   cutg->SetPoint(13,1.60242,16.6337);
   cutg->SetPoint(14,1.60304,12.3845);
   cutg->SetPoint(15,1.60304,7.41927);
   cutg->SetPoint(16,1.60398,2.40625);
   cutg->SetPoint(17,1.60367,-0.028646);
   cutg->SetPoint(18,1.52563,0.0668401);
   cutg->SetPoint(19,1.4523,0.0190971);
   cutg->SetPoint(20,1.44196,-0.0763891);
   cutg->SetPoint(21,1.44258,-0.171875);

   // taglio energia del fascio vs angolo theta del neutrone per selezionare gli eventi del canale di fotoproduzione di etaprimo su neutrone sopra la soglia
   TCutG *cutg2 = new TCutG("cutg2",29);
   cutg2->SetVarX("beam energy(GeV)");
   cutg2->SetVarY("Neutron theta(°)");
   cutg2->SetTitle("Graph");
   cutg2->SetFillStyle(1000);
   cutg2->SetPoint(0,1.56128,0.0469697);
   cutg2->SetPoint(1,1.56071,2.35352);
   cutg2->SetPoint(2,1.561,6.33755);
   cutg2->SetPoint(3,1.56214,10.2377);
   cutg2->SetPoint(4,1.56186,13.0056);
   cutg2->SetPoint(5,1.56157,16.1089);
   cutg2->SetPoint(6,1.56214,18.3316);
   cutg2->SetPoint(7,1.54982,17.7025);
   cutg2->SetPoint(8,1.53234,16.2347);
   cutg2->SetPoint(9,1.51401,14.6411);
   cutg2->SetPoint(10,1.4971,12.8378);
   cutg2->SetPoint(11,1.4822,10.9506);
   cutg2->SetPoint(12,1.47103,9.44091);
   cutg2->SetPoint(13,1.46014,7.1763);
   cutg2->SetPoint(14,1.44954,4.74394);
   cutg2->SetPoint(15,1.44438,2.8987);
   cutg2->SetPoint(16,1.44352,1.05346);
   cutg2->SetPoint(17,1.44323,-0.0788419);
   cutg2->SetPoint(18,1.47762,0.00503251);
   cutg2->SetPoint(19,1.53263,0.00503251);
   cutg2->SetPoint(20,1.56071,0.00503251);
   cutg2->SetPoint(21,1.56128,0.0469697);

   // taglio energia fascio vs angolo theta del protone per il canale di fotoproduzione di etaprimo su protone, allargato di un bin di energia sotto la soglia
   TCutG *cutg3 = new TCutG("cutg3",21);
   cutg3->SetVarX("beam energy(GeV)");
   cutg3->SetVarY("Proton Theta(°)");
   cutg3->SetTitle("Graph");
   cutg3->SetFillStyle(1000);
   cutg3->SetPoint(0,1.43309,-0.152027);
   cutg3->SetPoint(1,1.43309,1.97635);
   cutg3->SetPoint(2,1.43517,4.61149);
   cutg3->SetPoint(3,1.44378,7.5);
   cutg3->SetPoint(4,1.45417,9.72973);
   cutg3->SetPoint(5,1.46545,11.9088);
   cutg3->SetPoint(6,1.48118,14.3919);
   cutg3->SetPoint(7,1.49899,16.1655);
   cutg3->SetPoint(8,1.52036,18.0405);
   cutg3->SetPoint(9,1.54233,19.7635);
   cutg3->SetPoint(10,1.56697,21.2331);
   cutg3->SetPoint(11,1.58359,22.1959);
   cutg3->SetPoint(12,1.60051,22.5);
   cutg3->SetPoint(13,1.60081,20.3209);
   cutg3->SetPoint(14,1.6017,8.46284);
   cutg3->SetPoint(15,1.6011,0.050676);
   cutg3->SetPoint(16,1.5358,2.8657e-07);
   cutg3->SetPoint(17,1.47405,0.101352);
   cutg3->SetPoint(18,1.4325,-0.0506754);
   cutg3->SetPoint(19,1.48296,0.101352);
   cutg3->SetPoint(20,1.43309,-0.152027);

   // taglio energia fascio vs angolo theta del neutrone per il canale di fotoproduzione di etaprimo su neutrone, allargato di un bin di energia sotto la soglia
   TCutG *cutg4 = new TCutG("cutg4",20);
   cutg4->SetVarX("beam energy(GeV)");
   cutg4->SetVarY("Neutron theta(°)");
   cutg4->SetTitle("Graph");
   cutg4->SetFillStyle(1000);
   cutg4->SetPoint(0,1.56128,18.08);
   cutg4->SetPoint(1,1.561,14.3895);
   cutg4->SetPoint(2,1.56042,9.60866);
   cutg4->SetPoint(3,1.56014,4.74394);
   cutg4->SetPoint(4,1.56014,-0.0788419);
   cutg4->SetPoint(5,1.53607,-0.0369047);
   cutg4->SetPoint(6,1.49567,-0.0788419);
   cutg4->SetPoint(7,1.46673,0.0469697);
   cutg4->SetPoint(8,1.44352,-0.0788419);
   cutg4->SetPoint(9,1.43206,-0.0788419);
   cutg4->SetPoint(10,1.4332,1.55671);
   cutg4->SetPoint(11,1.4355,3.69551);
   cutg4->SetPoint(12,1.44323,5.91818);
   cutg4->SetPoint(13,1.44925,7.38598);
   cutg4->SetPoint(14,1.46157,9.48285);
   cutg4->SetPoint(15,1.47504,12.1668);
   cutg4->SetPoint(16,1.51687,16.6122);
   cutg4->SetPoint(17,1.54209,17.9961);
   cutg4->SetPoint(18,1.56128,18.2058);
   cutg4->SetPoint(19,1.56128,18.08);


   ofstream scrivimi("beamE.txt");

  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  for (Long64_t jentry=0; jentry<nentries; jentry++){
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    fproangles.clear();
    fdeuangles.clear();
    cphoton.clear();
    proton.clear();
    fneutron.clear();
    pionangles.clear();
    fphotonangles.clear();
    fDeuteron.clear();
    Tof_fpro.clear();
    Tof_fdeu.clear();

    fch=0;

    // if (Cut(ientry) < 0) continue;
    idrun = Idrun;
    idevt = Idevt;
    ipol  = int(Ipol);


    int fpro = 0;
    int cpro = 0;

    beam.SetPxPyPzE(0.,0.,Eg_tag_strip[0],Eg_tag_strip[0]);

    //fa un ciclo su tutte le traccie rivelate al centro
    for(int i=0; i<Nass_3; i++){
      if(Itipo_track[i]==11){ // le tracce che sono associate al codice 002(solo BGO) rappresentano eventi centrali neutri caratterizzati da una deposizione di energia di cluster minore 15 MeV(soglia software) ma comunque maggiore di 2 MeV su singolo cristallo(soglia hardware)
        TLorentzVector CandidatePhoton;
        CandidatePhoton.SetPxPyPzE(Eclusc_track[i]*sin(Thet_centr_track[i]/180.*M_PI)*cos(Phi_centr_track[i]/180.*M_PI),Eclusc_track[i]*sin(Thet_centr_track[i]/180.*M_PI)*sin(Phi_centr_track[i]/180.*M_PI),Eclusc_track[i]*cos(Thet_centr_track[i]/180.*M_PI),Eclusc_track[i]);
        cphoton.push_back(CandidatePhoton);
      }
      if(Itipo_track[i]==13 || Itipo_track[i]==14 || Itipo_track[i]==24)
      { // aggiungiamo le tracce classificate come Itipo_track=24 che sono eventi carichi associati al codice 112(camere a filo-barrel-BGO) con segnale sopra soglia sul barrel e tempo di volo all'interno della finestra temporale
        protondxE->Fill(Eclusc_track[i], Dedx_track[i]);
        int I_bar_track_ = static_cast<int>(I_bar_track[i]);
        protondxE_Ind_bar->Fill(I_bar_track_,Eclusc_track[i], Dedx_track[i]);
        TLorentzVector CandidateProton;
        if(ProtonCentrCut->IsInside(Eclusc_track[i], Dedx_track[i])){
          protondxEsoloPro->Fill(Eclusc_track[i], Dedx_track[i]);
          //   double Etotpro = Eclusc_track[i] + Dedx_track[i]/(2*sin(Thet_centr_track[i]/180.*M_PI)*1000.) + RMP;
          double Etotpro = Eclusc_track[i] + RMP;
          TLorentzVector CandidateProton;
          double Ppro = sqrt(Etotpro*Etotpro-RMP*RMP);
          CandidateProton.SetPxPyPzE(Ppro*sin(Thet_centr_track[i]/180.*M_PI)*cos(Phi_centr_track[i]/180.*M_PI), Ppro*sin(Thet_centr_track[i]/180.*M_PI)*sin(Phi_centr_track[i]/180.*M_PI), Ppro*cos(Thet_centr_track[i]/180.*M_PI),Etotpro);
          proton.push_back(CandidateProton);
          ccharged.push_back(CandidateProton);
          cpro++;
        }
        else if(PionCentrCut->IsInside(Eclusc_track[i], Dedx_track[i])){
          protondxEsoloPio->Fill(Eclusc_track[i], Dedx_track[i]);
          pair<double,double> tempangle;
          tempangle.first=Thet_centr_track[i];
          tempangle.second=Phi_centr_track[i];
          pionangles.push_back(tempangle);
        }
        // qui aggiungiamo l'informazione su un carico centrale che non è nè un pione nè un protone
        else{
           double Etotpro = Eclusc_track[i] + RMP;
          TLorentzVector CandidateProton;
          double Ppro = sqrt(Etotpro*Etotpro-RMP*RMP);
          CandidateProton.SetPxPyPzE(Ppro*sin(Thet_centr_track[i]/180.*M_PI)*cos(Phi_centr_track[i]/180.*M_PI), Ppro*sin(Thet_centr_track[i]/180.*M_PI)*sin(Phi_centr_track[i]/180.*M_PI), Ppro*cos(Thet_centr_track[i]/180.*M_PI),Etotpro);
          ccharged.push_back(CandidateProton);

        }
      }
      // identifichiamo i deuteroni al centro(taglio di prova: tutti gli eventi sopra la banana dei protoni sono nuclei di deuterio)
      /*
      if(cartella=="2005_d2"){
        if(FitDeuCentr_2005_d2->Eval(Eclusc_track[i])<=Dedx_track[i]){
           double Etotdeu = Eclusc_track[i] +RMD;
           double Pdeu = sqrt(Etotdeu*Etotdeu-RMD*RMD);
           TLorentzVector CandidateDeuteron;
           CandidateDeuteron.SetPxPyPzE(Pdeu*sin(Thet_centr_track[i]/180.*M_PI)*cos(Phi_centr_track[i]/180.*M_PI), Pdeu*sin(Thet_centr_track[i]/180.*M_PI)*sin(Phi_centr_track[i]/180.*M_PI), Pdeu*cos(Thet_centr_track[i]/180.*M_PI),Etotdeu);
           cDeuteron.push_back(CandidateDeuteron);
        }
      }
      */
    }
    // solo per verifica
    //cout<<"########"<<endl;
    /*
    for(int i=0; i<N_traf; i++){
      cout<<"Iass_trf["<<i<<"]="<<int(Iass_trf[i])<<endl;
    }
    */
    //ciclo su tutte le traccie delle particelle in avanti
    for(int i =0; i<Nparf; i++){
     //cout<<i<<"  "<<Nparf<<"----"<< int(Index_trf[Nparf])<<endl;
     //     cout<< int(Iass_trf[Index_trf[Nparf]])<<endl;//ci aspettiamo numeri pari a 1 - 7
     //
     //  cout <<int(Ind_traf[Index_trf[Nparf]][2])<<endl;

     N_show_->Fill(Nshow); // filliamo con il numero di colpi nel Russian Wall

     int index_trf_=0;
     if(Index_trf[i]!=0){
       index_trf_ = Index_trf[i]-1;// scaliamo il valore di Index_trf[i] di 1 perchè il Root parte da 0 e non da 1 come il Fortran
     }
      int index = int(Iass_trf[index_trf_]);

      ind_trf->Fill(index_trf_); // filliamo con l'indice assegnato alla traccia in base all'associazione dei rivelatori corrispondente
      Iass_trf_->Fill(index); // filliamo con l'indice asegnato alla traccia che identifica il tipo di particella


      // solo per verifica
      /*
      cout<<"Nparf="<<Nparf<<endl;
      cout<<"Index_trf["<<i<<"]="<<int(index_trf_)<<"  "<<endl;
      cout<<"Iass_trf["<<int(index_trf_)<<"]="<<int(Iass_trf[index_trf_])<<"  "<<endl;
      cout<<"index="<<index<<"  "<<endl;
      */
      //int indexRW = int(Ind_traf[Index_trf[i]][2]);
      if(index==1){//neutri
        //double DEDX_TRACCIAF   = Eshow[indexRW]*1000.*Cos_trf[i][2];
        //double TOF_F           = Tof_trf[i]*fabs(cos(Theta_trf[i])); // per normalizzare il tempo di volo bisogna moltiplicare per il coseno di theta della traccia
        double DEDX_F_N = De_trf[i]; // perdite di energia dei neutri in avanti(nel prean fanciullo vengono utilizzate le stesse varibili per neutri e carichi sulla peridta di energia e di tempo di volo)
        double TOF_F_N = Tof_trf[i]; // tempo di volo dei neutri in avanti
        ddxtof_f->Fill(TOF_F_N,DEDX_F_N);
        E_show_->Fill(Eshow[i]); // filliamo con la deposizione di energia nel russian wall
        Theta_trf_fn->Fill(Theta_trf[i]);
        Phi_trf_fn->Fill(Phi_trf[i]);
        if(Tof_trf[i]>=12){  //regno dei neutroni
          double beta     = 335./(Tof_trf[i]*CLIGHT*1.E-09);
          if((1-beta*beta)>0) {
            TLorentzVector CandidatefNeutron;
            double gamma = 1./sqrt(1-beta*beta);
            double ENE_FNeutron = gamma* RMN; //energia totale
            double Pfneu = sqrt(ENE_FNeutron*ENE_FNeutron-RMN*RMN);
            CandidatefNeutron.SetPxPyPzE(Pfneu*sin(Theta_trf[i]/180.*M_PI)*cos(Phi_trf[i]/180.*M_PI), Pfneu*sin(Theta_trf[i]/180.*M_PI)*sin(Phi_trf[i]/180.*M_PI), Pfneu*cos(Theta_trf[i]/180.*M_PI),ENE_FNeutron);
            fneutron.push_back(CandidatefNeutron);
          }
        }
        else if(Tof_trf[i]>=7.5 && Tof_trf[i]<=12.5){//regno fotoni ...solo gli angoli
          pair<double,double> tempangle;
          tempangle.first=Theta_trf[i];
          tempangle.second=Phi_trf[i];
          fphotonangles.push_back(tempangle);
        }
      }
      if(index!=1){//particelle cariche, quelle per le quali tutti i rivelatori sensibili ai carichi hanno sparato

        // if(index==6||index==7){//particelle cariche, quelle per le quali tutti i rivelatori sensibili ai carichi hanno sparato
        /*
        if(index==4){
          cout<<"non dovrei essere qui"<<endl;
        }
        */
        if(index==5){
                DIST_WALL = 335.; //in cm
        }
        else{
                // For tracks n.2 (0-1-0), 3 (0-1-1), 6 (1-1-0) and 7 (1-1-1),
                //the distance for tof is the OW
                DIST_WALL = 301.53;
        }
        Theta_trf_fch->Fill(Theta_trf[i]);
        Phi_trf_fch->Fill(Phi_trf[i]);
        ddxtof_fow->Fill(Tof_trf[i],De_trf[i]);
        if(cartella!="2005_d1"){  // nella cartella 2005_d1 non ci sono protoni rivelati in avanti
         if(ProtonForwCut->IsInside(Tof_trf[i], De_trf[i])){ //regno dei protoni in avanti
           pair<double,double> pangles;
           pangles.first = Theta_trf[i];
           pangles.second = Phi_trf[i];
           fproangles.push_back(pangles);
           Tof_fpro.push_back(Tof_trf[i]);
           double beta     = DIST_WALL/(Tof_trf[i]*CLIGHT*1.E-09);
	 if((1-beta*beta)>0){
	  TLorentzVector CandidatefProton;
	  double gamma = 1./sqrt(1-beta*beta);
	  double ENE_FPROTON = gamma* RMP; //energia totale
	  double Pfpro = sqrt(ENE_FPROTON*ENE_FPROTON-RMP*RMP);
	  CandidatefProton.SetPxPyPzE(Pfpro*sin(Theta_trf[i]/180.*M_PI)*cos(Phi_trf[i]/180.*M_PI), Pfpro*sin(Theta_trf[i]/180.*M_PI)*sin(Phi_trf[i]/180.*M_PI), Pfpro*cos(Theta_trf[i]/180.*M_PI),ENE_FPROTON);
	  proton.push_back(CandidatefProton);
	  fpro++;
	 }
         }
         if(PionForwCut->IsInside(Tof_trf[i], De_trf[i])){ //regno dei pioni in avanti
	  pair<double,double> tempangle;
	  tempangle.first=Theta_trf[i];
	  tempangle.second=Phi_trf[i];
	  pionangles.push_back(tempangle);
         }
         /*
         if(cartella!="2005_d2"){
           if(De_trf[i]>=0.007 && Tof_trf[i]>10.44){
             Tof_fch.push_back(Tof_trf[i]);
             fch++;
           }
         }
         if(cartella!="1999_d2"){
           if(De_trf[i]>=0.008 && Tof_trf[i]>=11){
             Tof_fch.push_back(Tof_trf[i]);
             fch++;
           }
         }
         if(cartella!="1999_d1"){
           if(De_trf[i]>=0.012 && Tof_trf[i]<=12.45){
             Tof_fch.push_back(Tof_trf[i]);
             fch++;
           }
         }
         if(cartella!="2001_d"){
           if(De_trf[i]>=0.007 && Tof_trf[i]<=10.84){
             Tof_fch.push_back(Tof_trf[i]);
             fch++;
           }
         }
         if(cartella!="2006_d"){
           if(De_trf[i]>=0.007 && Tof_trf[i]<=13.05){
             Tof_fch.push_back(Tof_trf[i]);
             fch++;
           }
         }
         if(cartella!="2002_d1"){
           if(De_trf[i]>=0.009 && Tof_trf[i]<=12.34){
             Tof_fch.push_back(Tof_trf[i]);
             fch++;
           }
         }
         if(cartella!="2002_d3"){
           if(De_trf[i]>=0.010 && Tof_trf[i]<=11.14){
             Tof_fch.push_back(Tof_trf[i]);
             fch++;
           }
         }
         if(cartella!="2002_d2"){
           if(De_trf[i]>=0.005 && Tof_trf[i]<=9.76){
             Tof_fch.push_back(Tof_trf[i]);
             fch++;
           }
         }
         */
         if(DeuForwCut!=0){
           if(DeuForwCut->IsInside(Tof_trf[i], De_trf[i])){ // regno dei deutoni in avanti
	    double beta = DIST_WALL/(Tof_trf[i]*CLIGHT*1.E-09);
	    pair<double,double> dangles;
	    dangles.first= Theta_trf[i];
	    dangles.second = Phi_trf[i];
	    fdeuangles.push_back(dangles);
	    Tof_fdeu.push_back(Tof_trf[i]);
	    if((1-beta*beta)>0){
	     TLorentzVector CandidatefDeu;
	     double gamma = 1./sqrt(1-beta*beta);
	     double ENE_FDEUTERON = gamma* RMD; // energia totale
	     double Pfdeu = sqrt(ENE_FDEUTERON*ENE_FDEUTERON-RMD*RMD);
	     //double Pfdeu = sqrt(ENE_FDEUTERON*ENE_FDEUTERON-RMP*RMP);
	     CandidatefDeu.SetPxPyPzE(Pfdeu*sin(Theta_trf[i]/180.*M_PI)*cos(Phi_trf[i]/180.*M_PI), Pfdeu*sin(Theta_trf[i]/180.*M_PI)*sin(Phi_trf[i]/180.*M_PI),Pfdeu*cos(Theta_trf[i]/180.*M_PI),ENE_FDEUTERON);
	     fDeuteron.push_back(CandidatefDeu);
	    }
            //Tof_fch.push_back(Tof_trf[i]);
            //fch++;
            }
         }
         /*
         if(!PionForwCut->IsInside(Tof_trf[i], De_trf[i]) && !ProtonForwCut->IsInside(Tof_trf[i], De_trf[i]) && !DeuForwCut->IsInside(Tof_trf[i], De_trf[i])){
           //Tof_fch.push_back(Tof_trf[i]);
           //fch++;
         }
         */
       }
      }
     }
    //tree->Fill();
    /*
    // taglio Eg vs Theta protone per filtrare gli eventi del canale di fotoproduzione di etaprimo su protone(e ridurre il peso delle rootople)
    if(proton.size()==1 && fneutron.size()==0 && fDeuteron.size()==0){
     if(beam.E()>=1.447){
      TLorentzVector *p=proton.data();
      if(cutg->IsInside(beam.E(),p[0].Theta()*TMath::RadToDeg())){
       scrivimi<<beam.E()<<endl;
       tree->Fill();
      }
     }
    }
    */
      /*
      if(cutg->IsInside(beam.E(),p[0].Theta()*TMath::RadToDeg())){
        tree->Fill();
      }
      */
    // taglio che seleziona soltanto gli eventi con un protone e 4 fotoni centrali
    /*
    if(proton.size()==1 && cphoton.size()==4){
     TLorentzVector *p=proton.data();
     tree->Fill();
    }
    */
    /*
    if(fneutron.size()==1 && cphoton.size()>=2){
      scrivimi<<beam.E()<<endl;
      tree->Fill();
    }
    */
    /*
     // taglio Eg vs Theta neutrone per filtrare gli eventi del canale di fotoproduzione di etaprimo su neutrone legato nel deuterio(per ridurre il peso delle rootople)
    if(fneutron.size()==1){
      TLorentzVector *n=fneutron.data();
      if(cutg2->IsInside(beam.E(),n[0].Theta()*TMath::RadToDeg())){
        tree->Fill();
      }
    }
    */
     /*
    // taglio che seleziona soltanto gli eventi con un deutone come unica particella carica emessa in avanti(la rootopla in output verrà analizzato per investigare i canali di fotoproduzione coherent di Omega su deuterio e rho su deuterio)
    if(fDeuteron.size()==1 && proton.size()==0 && fneutron.size()==0){
      tree->Fill();
    }
    */

    /* // selezione dei soli eventi con un deutone in avanti nello stato finale  e 2 fotoni centrali(la rootopla in output servirà ad investigare il canale di fotoproduzione coherente di eta su deuterio)
    if(fDeuteron.size()==1 && proton.size()==0 && fneutron.size()==0 && cphoton.size()==2){
      tree->Fill();
    }
    */

    // selezione dei soli eventi con un deutone i avanti nello stato finale e 4 fotoni centrali(la rootopla in output verrà analizzata per investigare i canali di fotoproduzione di eta e pi0 su deuterio e doppio pi0 su deuterio)
    if(fDeuteron.size()==1 && proton.size()==0 && fneutron.size()==0 && cphoton.size()==4){
      tree->Fill();
    }


    // taglio che seleziona soltanto gli eventi con un protone in avanti e 4 fotoni centrali con la  sostituzione della massa del protone rivelato con quella del deuterio
    /*
    if(proton.size()==1 && fneutron.size()==0 && cphoton.size()==4 && fDeuteron.size()==0){
      tree->Fill();
    }
    */
    // selezioniamo soltanto gli eventi con una particella carica in avanti e 4 particelle neutre centrali
    /*
    if(fch==1 && cphoton.size()==4 && fneutron.size()==0 && fphotonangles.size()==0 ){
      tree->Fill();
    }
    */
    //if(cphoton.size()==4 && pionangles.size()==0 && fneutron.size()==0 && ((fch+ccharged.size())==1)){
      //tree->Fill();
    //}
    /*
    // selezione eventi di con 4 fotoni centrali e un protone e un neutrone in avanti(reazione di break-up du deuterio con produzione di 2 pioni neutri)
    if(cphoton.size()==4 && fneutron.size()==1 && proton.size()==1){
      tree->Fill();
    }
    */
    /*
    if(cphoton.size()==4 && fDeuteron.size()==1){
       tree->Fill();
    }
    */

     //tree->Fill(); // fotografa il tree
     //scrivimi<<beam.E()<<endl;
  }

  scrivimi.close();

  // i due blocchi successivi pervono  a controllare i tagli che identificano i protoni e  i pioni in avanti e centrali
  TCanvas *c1 = new TCanvas("c1","",700,600);
  ProtonCentrCut->SetLineColor(2);
  PionCentrCut->SetLineColor(1);
  protondxE->SetName("De_dx vs Er");
  protondxE->Draw("colz");
  protondxE->GetXaxis()->SetTitleSize(0.063);
  protondxE->GetYaxis()->SetTitleSize(0.063);
  protondxE->GetXaxis()->CenterTitle(1);
  protondxE->GetXaxis()->SetTitleOffset(0.7);
  protondxE->GetYaxis()->SetTitleOffset(0.7);
  protondxE->GetXaxis()->SetTitle("Er(GeV)");
  protondxE->GetYaxis()->SetTitle("Dedx(GeV/cm)");
  ProtonCentrCut->SetLineColor(2);
  PionCentrCut->SetLineColor(5);
  ProtonCentrCut->Draw("same");
  PionCentrCut->Draw("same");

  TLegend *l = new TLegend(0.5,0.92,0.9,1.0);
  l->AddEntry(ProtonCentrCut,"Protons","l");
  l->AddEntry(PionCentrCut,"Pions","l");
  l->Draw("same");
  string sep="/home/thinkpadinfn/GRAAL/";
  string Name1 ="De_dx_vs_Er_"+cartella+ ".pdf";
  c1->SaveAs(Name1.c_str());

  if(DeuForwCut!=0){
  TCanvas *c3 = new TCanvas("c3",str3.c_str(),700,600);
  c3->SetLogz();
  ProtonForwCut->SetLineColor(1);
  PionForwCut->SetLineColor(5);
  DeuForwCut->SetLineColor(2);
  ProtonForwCut->SetLineWidth(5);
  PionForwCut->SetLineWidth(5);
  DeuForwCut->SetLineWidth(5);
  ddxtof_fow->Draw("colz");
  ddxtof_fow->SetName("De_dx vs Tof");
  ddxtof_fow->GetXaxis()->SetTitle("Tof(ns)");
  ddxtof_fow->GetYaxis()->SetTitle("Dedx(MeV/cm)");
  ddxtof_fow->GetXaxis()->SetTitleSize(0.063);
  ddxtof_fow->GetYaxis()->SetTitleSize(0.063);
  ddxtof_fow->GetXaxis()->CenterTitle(1);
  ddxtof_fow->GetYaxis()->CenterTitle(1);
  ddxtof_fow->GetXaxis()->SetTitleOffset(0.9);
  ddxtof_fow->GetYaxis()->SetTitleOffset(0.9);
  PionForwCut->Draw("same");
  ProtonForwCut->Draw("same");
  DeuForwCut->Draw("same");


  TLegend *l2 = new TLegend(0.5,0.76,0.9,0.86);
  l2->AddEntry(ProtonForwCut,"Protons","l");
  l2->AddEntry(PionForwCut,"Pions","l");
  l2->AddEntry(DeuForwCut,"Deutons","l");
  l2->Draw("same");
  string Name2 = "De_dx_vs_Tof_charged_"+cartella+".pdf";
  c3->SaveAs(Name2.c_str());
  }

  TCanvas *c4 = new TCanvas("c4","",700,600);
  ddxtof_f->GetXaxis()->CenterTitle();
  ddxtof_f->GetXaxis()->SetTitleSize(0.063);
  ddxtof_f->GetXaxis()->SetTitle("Tof(ns)");
  ddxtof_f->GetYaxis()->SetTitleSize(0.063);
  ddxtof_f->GetYaxis()->SetTitle("De/dx(MeV/cm)");
  ddxtof_f->GetYaxis()->CenterTitle();
  ddxtof_f->Draw("");
  string Name3 ="De_dx_vs_Tof_neutr"+cartella+".pdf";
  c4->SaveAs(Name3.c_str());


  /*
  TFile *fileI = new TFile("../Tracks_spectra.root","update");
  protondxE->Write();
  ddxtof_f->Write();
  ddxtof_fow->Write();
  fileI->Close();
  */

  /*
  TFile *fileII = new TFile("../Tracks_spectra_forward.root","update");
  ddxtof_fow_tracks_5_6->Write();
  ddxtof_fow_tracks_2_3->Write();
  ddxtof_fow_tracks_4->Write();
  fileII->Close();
  */

  /*
  TFile *fileIII = new TFile("../Centr_Tracks_dedx_Er_Ind_bar_1998_uv.root","recreate");
  protondxE_Ind_bar->Write();
  fileIII->Close();
  */

  /*
  TFile *fileIV = new TFile("Forward_photons_1998_uv.root","recreate");
  ddxtof_f->Write();
  fileIV->Close();
  */

  /*
  TFile *fileV = new TFile(str12.c_str(),"recreate");
  ddxtof_fow->Write();
  Theta_trf_fch->Write();
  Phi_trf_fch->Write();
  N_show_->Write();
  E_show_->Write();
  Iass_trf_->Write();
  ind_trf->Write();
  ddxtof_f->Write();
  Theta_trf_fn->Write();
  Phi_trf_fn->Write();

  fileV->Close();
  */

  dati->Write();
  dati->Close();


}
   



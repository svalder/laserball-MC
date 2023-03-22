// ---------------------------------------------------------
// Goal:          Simulate photon tracking through laserball diffuser flask
// Authors:        Sammy Valder and Martti Nirkko, (2022)
// Compile & run: clear && g++ -g -o diffuser.exe diffuser.C `root-config --cflags --libs` && ./diffuser.exe
// ---------------------------------------------------------

// Helper functions (includes everything else)
#include "include/utilities.C"
#include "include/config.C"

// *****************************************************************************
// Function declarations
TVector3 diffuser(double& tracklen, double scatlen);
double distance_to_wall(TVector3&, TVector3&);
double reflection_prob(const double& n1, const double& n2, const double &impact);
TVector3 propagate(TVector3&, TVector3&, double& distance);
TVector3 scatter(TVector3&, TVector3&, double& impact);
TVector3 reflect_or_refract(TVector3&, TVector3&, double& impact, bool& exitflask);
double GetScatteringLength(double density);

// Other global objects (TODO - BAD PRACTICE)
TRandom3* gen = new TRandom3();
TF1* aperture = new TF1("aperture","cos(pi/2*x/[0])",0,NA);
TVector3 pos, dir, newpos, newdir, endpos, enddir;
TGraph trackS, trackT, trackU; // side view, top view
TH1D htrk("htrk","Single photon tracking (statistics);Distance between scatters [mm];Events",NBINS/2,0,5);
TFile outfile("diffuser.root","RECREATE");
int event;

// *****************************************************************************
// Main program
int main(int argc, char** argv) {

  // Set output verbosity (kPrint, kInfo, kWarning, kError, kBreak, kSysError, kFatal)
  gErrorIgnoreLevel = kWarning;

  // Scattering length of photons in diffuser medium
  double scatlen;

  // Initialisation
  gen->SetSeed(0);
  aperture->SetParameter(0,NA);
  TVector3 endpos;
  TVector3 outdir;
  double length;
  gStyle->SetTitleOffset(1.1,"x");
  gStyle->SetTitleOffset(1.4,"y");
  std::string ext;
  std::string tlen;
  std::string tphi;  
  std::string tcth;
  std::string tang;
  std::string ttpo;
  std::string ttag;

  // Loop over all offsets
  for (int it=0; it<sizeof(var)/sizeof(double); it++) {
    
    if (VAR == 0){
      OFFSET = var[it];
      tlen = Form("hesc_z0_%.0fmm",OFFSET*1e3);
      tphi = Form("hphi_z0_%.0fmm",OFFSET*1e3);
      tcth = Form("hcth_z0_%.0fmm",OFFSET*1e3);
      tang = Form("hang_z0_%.0fmm",OFFSET*1e3);
      ttpo = Form("htpo_z0_%.0fmm",OFFSET*1e3);
      ttag = Form("htag_z0_%.0fmm",OFFSET*1e3);
      cout << endl << "SIMULATING DIFFUSER WITH INJECTION POINT OFFSET Z = " << OFFSET << " mm" << endl;
    }

    if (VAR == 1){
      con_bub = var[it];
      tlen = Form("hesc_z0_%.0fmgmL",con_bub*1e6);
      tphi = Form("hphi_z0_%.0fmgmL",con_bub*1e6);
      tcth = Form("hcth_z0_%.0fmgmL",con_bub*1e6);
      tang = Form("hang_z0_%.0fmgmL",con_bub*1e6);
      ttpo = Form("htpo_z0_%.0fmgmL",con_bub*1e6);
      ttag = Form("htag_z0_%.0fmgmL",con_bub*1e6);
      cout << endl << "SIMULATING DIFFUSER WITH GLASS BEAD CONCENTRATION = " << con_bub*1e3 << " mg/ml" << endl;
    }

    scatlen = 10*GetScatteringLength(con_bub); // mean free path [mm]
    if(VERBOSE > 1) cout << "Scattering length = " << scatlen << " mm" << endl;

    TH1D hesc(tlen.c_str(),"Escape time from diffuser;t [ns];Events [#times 10^{3}]",10*NBINS,0,10.);
    TH1D hphi(tphi.c_str(),Form("Azimuthal distribution;#phi [#pi];Events [#times 10^{3}]"),NBINS,-1,1);
    TH1D hcth(tcth.c_str(),Form("Polar distribution;cos(#theta) [ ];Events [#times 10^{3}]"),NBINS,-1,1);
    TH2D hang(tang.c_str(),"Angular distribution of diffuser;#phi [#pi];cos(#theta) [ ]",NBINS/2,-1,1,NBINS/2,-1,1);
    TH2D htpo(ttpo.c_str(),"Timing as a function of polar angle", NBINS,-1,1,NBINS,0,3);
    TH3D htag(ttag.c_str(),"Temporal distribution across diffuser;#phi [#pi];cos(#theta);Escape time [ns]",NBINS/2,-1,1,NBINS/2,-1,1,NBINS/2,0,10);
    
    // Generate events
    for (event=0; event<NEVENTS; event++) {
      printProgress(event,NEVENTS);
      endpos = diffuser(length, scatlen);
      outdir = endpos.Unit();
      hesc.Fill(length/cs); // assume constant speed in silicone, neglect time spend in bubbles
      hphi.Fill(outdir.Phi()/pi);
      hcth.Fill(cos(outdir.Theta()));
      hang.Fill(outdir.Phi()/pi,cos(outdir.Theta()));
      htpo.Fill((cos(outdir.Theta())),length/cs);
      htag.Fill(outdir.Phi()/pi,cos(outdir.Theta()),length/cs);
    }

    // Plot distributions
    double scale = 1e-3;
    TCanvas c("c","",1200,1200);
    c.Divide(2,2);
    c.cd(1)->SetGrid();
    hphi.Scale(scale);
    hphi.SetMinimum(0);
    //hphi.SetAxisRange(0,1.25*NEVENTS/NBINS*scale,"Y");
    //hphi.GetYaxis()->SetTitleOffset(1.2);
    hphi.Draw();
    c.cd(2)->SetGrid();
    hcth.Scale(scale);
    hcth.SetMinimum(0);
    //hcth.SetAxisRange(0,1.25*NEVENTS/NBINS*scale,"Y");
    //hcth.GetYaxis()->SetTitleOffset(1.2);
    hcth.Draw();
    c.cd(3)->SetGrid();
    hang.Draw("colz");
    hang.SetMinimum(0);
    //hang.SetAxisRange(0,5*NEVENTS/NBINS/NBINS,"Z");
    //hang.GetYaxis()->SetTitleOffset(1.2);
    hang.Draw("colz same");
    c.cd(4)->SetGrid();
    hesc.Scale(scale);
    hesc.SetMinimum(0);
    //hesc.SetAxisRange(0,5*NEVENTS/NBINS*scale,"Y");
    hesc.Draw();
    if (VAR==0){
      c.Print(Form("diffuser_z0_%.0fmm.png",OFFSET*1e3));
      c.Print(Form("diffuser_z0_%.0fmm.pdf",OFFSET*1e3));
    }
    if (VAR==1){
      c.Print(Form("diffuser_z0_%.2fmgmL.png",con_bub*1e3));
      c.Print(Form("diffuser_z0_%.2fmgmL.pdf",con_bub*1e3));
    }
    c.Close();

    TCanvas c2("c2","",1200,1200);
    htpo.Draw("colz");
    if (VAR==0){
      c2.Print(Form("diffuser_z0_timimngPolar_%.0fmm.png",OFFSET*1e3));
      c2.Print(Form("diffuser_z0_timimngPolar_%.0fmm.pdf",OFFSET*1e3));
    }
    if (VAR==1){
      c2.Print(Form("diffuser_z0_timimngPolar_%.2fmgmL.png",con_bub*1e3));
      c2.Print(Form("diffuser_z0_timimngPolar_%.2fmgmL.pdf",con_bub*1e3));
    }
    c2.Close();

    TCanvas *c3 = new TCanvas();
    htag.Draw("ISO");
    if (VAR==0){
      c3->Print(Form("tempAngDist_%.0fmm.png",OFFSET*1e3));
      c3->Print(Form("tempAngDist_%.0fmm.pdf",OFFSET*1e3));
    }
    if (VAR==1){
      c3->Print(Form("tempAngDist_%.2fmgmL.png",con_bub*1e3));
      c3->Print(Form("tempAngDist_%.2fmgmL.pdf",con_bub*1e3));
    }
    c3->Close();

    cout << "Time spread is " << getFWHM(&hesc) << " ns (FWHM)." << endl;
    
    outfile.Write();
  }
  outfile.Close();

  return 0;
}

// *****************************************************************************
// Track a single photon through the diffuser flask
// Inputs: none
// Output: escape direction of photon, seen from 6m sphere
TVector3 diffuser(double& tracklen, double scatlen) {
  
  // Photon position at injection point
  pos = e1;
  double A = pi*pow(roddiam/2.,2)*gen->Rndm(); // random within area
  pos.SetMag(sqrt(A/pi)); // radius [mm]
  pos.SetPhi(2*pi*gen->Rndm()); // angle [rad]
  pos += OFFSET*e3; // offset in z-direction

  // Photon direction at injection point
  dir = e3;
  dir.SetPhi(2*pi*gen->Rndm()); // azimuthal direction [0,2pi)
  double theta = aperture->GetRandom();
  dir.SetTheta(pi-theta); // downwards zenith direction [0,NA)
  
  // Photon tracking
  int step = 0;
  tracklen = 0;
  trackS.Set(0);
  trackT.Set(0);
  trackU.Set(0);
  // Injection point
  if (!event) {
    trackS.SetPoint(step,pos.X(),pos.Z());
    trackT.SetPoint(step,pos.X(),pos.Y());
    trackU.SetPoint(step,pos.Y(),pos.Z());
  }
  if (VERBOSE > 1) printf("ENTERING diffuser at (%8.3f %8.3f %8.3f), R=%6.3f\n",pos.X(),pos.Y(),pos.Z(),pos.Mag());
  bool exitflask = false;
  while (!exitflask) {
    if (VERBOSE > 2) printf("position (%8.3f %8.3f %8.3f), direction = (%8.3f %8.3f %8.3f)\n",pos.X(),pos.Y(),pos.Z(),dir.X(),dir.Y(),dir.Z());
    double y = gen->Rndm();
    double dsel = -scatlen*log(y);  // randomly selected distance in diffuser [mm]
    double dwall = distance_to_wall(pos,dir); // distance to flask wall [mm]
    double impact = sqrt(y);
    if (dsel < dwall) {
      // propagate selected distance
      newpos = propagate(pos,dir,dsel);
      tracklen += dsel;
      if (!event) htrk.Fill(dsel);
      // scatter in diffuser
      newdir = scatter(newpos,dir,impact);
      step++;
    } else {
      // propagate up to wall
      newpos = propagate(pos,dir,dwall);
      tracklen += dwall;
      if (!event) htrk.Fill(dwall);
      // exit or reflect internally
      newdir = reflect_or_refract(newpos,dir,impact,exitflask);
    }
    pos = newpos;
    dir = newdir;
    if(pos.Mag()>R+1e-6) printf("*** WARNING *** Position is (%8.3f %8.3f %8.3f), R=%6.3f\n",pos.X(),pos.Y(),pos.Z(),pos.Mag());
    if(fabs(dir.Mag()-1.)>1e-6) printf("*** WARNING *** Direction is (%8.3f %8.3f %8.3f), R=%6.3f\n",dir.X(),dir.Y(),dir.Z(),dir.Mag());
    if (!event) {
      trackS.SetPoint(step,pos.X(),pos.Z());
      trackT.SetPoint(step,pos.X(),pos.Y());
      trackU.SetPoint(step,pos.Y(),pos.Z());
    }
  }
  
  if (VERBOSE > 1) printf("EXITING diffuser after %d scatters at (%8.3f %8.3f %8.3f), R=%6.3f\n",step,pos.X(),pos.Y(),pos.Z(),pos.Mag());
  double RMAX = 6e3; // AV radius [mm] //TODO: Can change this for PMT settings
  endpos = propagate(pos,dir,RMAX);
  if (VERBOSE > 1) printf("End point (%8.3f %8.3f %8.3f), R=%6.3f\n",endpos.X(),endpos.Y(),endpos.Z(),endpos.Mag());
  if (!event) {
    trackS.SetPoint(step+1,endpos.X(),endpos.Z());
    trackT.SetPoint(step+1,endpos.X(),endpos.Y());
    trackU.SetPoint(step+1,endpos.Y(),endpos.Z());
    
    double LIM=75;
    gStyle->SetOptStat(0);
    //gStyle->SetAxisLabelOffset?
    TCanvas c("c","",1200,1200);
    c.Divide(2,2);
    c.cd(1)->DrawFrame(-LIM,-LIM,LIM,LIM,"Single photon tracking (front view);X [mm];Z [mm]");
    // Flask
    TEllipse ball(0,0,R);
    ball.Draw("L same");
    // Rod
    TBox rod(-roddiam/2,OFFSET,roddiam/2,LIM); // must be within limits of frame
    rod.SetFillColor(0);
    rod.SetLineColor(4);
    rod.Draw("L same");
    // Photon track
    trackS.SetLineColor(2);
    trackS.Draw("L same");
    c.cd(3)->DrawFrame(-LIM,-LIM,LIM,LIM,"Single photon tracking (top view);X [mm];Y [mm]");
    // Flask
    ball.Draw("L same");
    // Rod
    TEllipse rodc(0,0,roddiam/2);
    rodc.SetFillColor(0);
    rodc.SetLineColor(4);
    rodc.Draw("L same");
    // Photon track
    trackT.SetLineColor(2);
    trackT.Draw("L same");
    c.cd(2)->DrawFrame(-LIM,-LIM,LIM,LIM,"Single photon tracking (side view);Y [mm];Z [mm]");
    // Flask
    ball.Draw("L same");
    // Rod
    rod.Draw("L same");
    // Photon track
    trackU.SetLineColor(2);
    trackU.Draw("L same");
    c.cd(4);
    htrk.GetXaxis()->SetTitleOffset(1.1);
    htrk.GetYaxis()->SetTitleOffset(1.4);
    htrk.Draw();
    c.Print(Form("photon_track_z0=%.1f.png",OFFSET));
    c.Print(Form("photon_track_z0=%.1f.pdf",OFFSET));
    c.Close();
  }

  return endpos;
}

// *****************************************************************************
// Calculate distance to inner flask boundary
double distance_to_wall(TVector3& pos, TVector3& dir) {
  double x = pos.X();
  double y = pos.Y();
  double z = pos.Z();
  double alpha = dir.X();
  double beta  = dir.Y();
  double gamma = dir.Z();
  double dist = -(x*alpha+y*beta+z*gamma)+sqrt(pow(x*alpha+y*beta+z*gamma,2)-(x*x+y*y+z*z-R*R));
  return dist;
}

// *****************************************************************************
// Reflection probability
double reflection_prob(const double& n1, const double& n2, const double& b) {
    // Consider transition between materials (n2 > n1)
    double ni = n1/n2;
    // Calculate Fresnel coefficients
    double rs = (sqrt(1.-b*b)-sqrt(ni*ni-b*b))/(sqrt(1.-b*b)+sqrt(ni*ni-b*b));
    double rp = (sqrt(ni*ni-b*b)-ni*ni*sqrt(1.-b*b))/(sqrt(ni*ni-b*b)+ni*ni*sqrt(1.-b*b));
    // Average over Fresnel coefficients:
    return 0.5*(rs*rs+rp*rp);
}

// *****************************************************************************
// Propagate through silicone
TVector3 propagate(TVector3& pos, TVector3& dir, double& dist) {
  double x = pos.X()+dir.X()*dist;
  double y = pos.Y()+dir.Y()*dist;
  double z = pos.Z()+dir.Z()*dist;
  return TVector3(x,y,z);
}

// *****************************************************************************
// Scatter off a glass bubble
TVector3 scatter(TVector3& pos, TVector3& dir, double& impact) {
  
  // Directional cosines
  double alpha = dir.X();
  double beta  = dir.Y();
  double gamma = dir.Z();
  
  // Impact parameter
  double b = impact;
  
  // Forward scattering angle depends on scattering modes A, B, C
  double costh;
  if (ng*b > 1) {
    // (A) forward scattering on glass bubble
    costh = 2.*b*b-1.;
  } else {
    // (B) scattering through glass bubble cavity
    costh = 2.*pow(ns*b*b+sqrt((1.-b*b)*(1.-ns*ns*b*b)),2)-1.;
    // (C) probability of reflecting despite nb <= 1
    double prob = reflection_prob(na,ng,b);
    double roll = gen->Rndm();
    if (roll < prob) costh = 2.*b*b-1.;
  }
  
  // Azimuthal angle chosen randomly [0,2*pi)
  double phi = 2.*pi*gen->Rndm();
  
  // New direction transformed back into global frame (R. Ford's thesis)
  double param = sqrt((1.-costh*costh)/(1.-gamma*gamma));
  double newalpha = alpha*costh + param*(alpha*gamma*cos(phi)-beta*sin(phi));
  double newbeta  = beta*costh + param*(beta*gamma*cos(phi)+alpha*sin(phi));
  double newgamma = gamma*costh - param*(1.-gamma*gamma)*cos(phi);  
  TVector3 scatdir(newalpha,newbeta,newgamma);
  return scatdir.Unit();
}

// *****************************************************************************
// Reflect at flask wall, or refract out of diffuser
TVector3 reflect_or_refract(TVector3& pos, TVector3& dir, double& impact, bool& exitflask) {

  // Photon position and directional cosines
  double x = pos.X();
  double y = pos.Y();
  double z = pos.Z();
  double alpha = dir.X();
  double beta  = dir.Y();
  double gamma = dir.Z();
  
  // Forward scattering angle depends on scattering modes A, B, C
  double prod = dir.Unit().Dot(pos.Unit());
  double sinOmega = ng*ng*(1.-fabs(prod));
  TVector3 scatdir;
  if (sinOmega > 1) {
    // (A) internal reflection
    scatdir = dir-2.*prod*pos;
    exitflask = false;
  } else {
    // (B) probability of reflecting despite sinOmega <= 1
    double prob = reflection_prob(ns,ng,impact);
    double roll = gen->Rndm();
    if (roll < prob) {
      scatdir = dir-2.*prod*pos;
      exitflask = false;
      if (VERBOSE > 1) printf("REFLECTING on diffuser wall at (%8.3f %8.3f %8.3f), R=%6.3f\n",x,y,z,pos.Mag());
    }
    // (C) refract out of the diffuser
    else {
      
      // Change of plan: Find rotation angles for simple frame
      double rot_X, rot_Z;
      GetRotationAngles(pos.Unit(),rot_Z,rot_X);
      
      // Rotate original position to simple frame
      TVector3 rotpos = pos.Unit();
      rotpos.RotateZ(-rot_X);
      rotpos.RotateY(-rot_Z);
      
      // Rotate original direction to simple frame
      TVector3 rotdir = dir;
      rotdir.RotateZ(-rot_X);
      rotdir.RotateY(-rot_Z);
      
      // Rotate so that everything is in x-z plane (simple frame 2D)
      double angZ = atan(rotdir.Y()/rotdir.X());
      rotdir.RotateZ(-angZ);
      
      // Calculate direction after refraction
      double u = pos.Unit()*dir.Unit(); // cos(theta) in both frames)
      double alp = sqrt(1-pow(ns/ng,2)*(1-u*u));
      double bet = ns/ng*sqrt(1-u*u);
      if (bet*rotdir.X() < 0) bet*=-1; // change sign
      scatdir.SetXYZ(bet,0,alp);
      
      // Check that rotated vectors are all in x-z plane
      if (VERBOSE > 2) {
        printf("-----\n");
        printf("rotpos = (%8.3f %8.3f %8.3f), length=%8.3f\n",rotpos.X(),rotpos.Y(),rotpos.Z(),rotpos.Mag());
        printf("rotdir = (%8.3f %8.3f %8.3f), length=%8.3f\n",rotdir.X(),rotdir.Y(),rotdir.Z(),rotdir.Mag());
        printf("newdir = (%8.3f %8.3f %8.3f), length=%8.3f\n",scatdir.X(),scatdir.Y(),scatdir.Z(),scatdir.Mag());
      }
      
      // Rotate back last step
      rotdir.RotateZ(angZ);
      scatdir.RotateZ(angZ);
      
      // Rotate back to original frame
      rotpos.RotateUz(pos.Unit());
      rotdir.RotateUz(pos.Unit());
      scatdir.RotateUz(pos.Unit());
      
      // Check output vectors
      if (VERBOSE > 2) {
        printf("-----\n");
        printf("rotpos = (%8.3f %8.3f %8.3f), length=%8.3f\n",rotpos.X(),rotpos.Y(),rotpos.Z(),rotpos.Mag());
        printf("rotdir = (%8.3f %8.3f %8.3f), length=%8.3f\n",rotdir.X(),rotdir.Y(),rotdir.Z(),rotdir.Mag());
        printf("newdir = (%8.3f %8.3f %8.3f), length=%8.3f\n",scatdir.X(),scatdir.Y(),scatdir.Z(),scatdir.Mag());
        printf("-----\n");
      }
      
      // Verify output vectors are correct
      if (VERBOSE > 2) {
        printf("CHECK RESULTS\n");
        double d = sqrt(1.-pow(ns/ng,2)*(1.-pow(fabs(prod),2)));
        printf("(1) Snell's law: 0 = %e\n",scatdir*pos.Unit()-d);
        printf("(2) Unit length: 0 = %e\n",scatdir.Mag()-1);
        printf("(3) All planar:  0 = %e\n",scatdir.Dot(pos.Cross(dir)));
        printf("-----\n");
      }
      
      // Set break condition
      exitflask = true;
    }
  }
  
  return scatdir.Unit();
}

// *****************************************************************************
// Get scattering length from theory, given a bubble mass density
double GetScatteringLength(double density) {
  if (VERBOSE > 0) cout << "Provided mass density of glass bubbles: " << density << " g/cm^3" << endl;
  
  // Scattering cross section of glass microspheres
  // https://en.wikipedia.org/wiki/Anomalous_diffraction_theory
  // https://www.gfdl.noaa.gov/wp-content/uploads/files/user_files/pag/lecture2008/lecture3.pdf
  //double n = ns/na; // refractive index (outer vs inner)
  //double p = 4.*pi*r_bub*(n-1.)/lambda; // phase delay of wave through sphere
  //double Q = 2. - 4./p*sin(p) + 4./(p*p)*(1.-cos(p)); // scattering efficiency (same as extinction efficiency)
  double Q = 1; // artificial (ignore ADT)
  double sigma = pi*r_bub*r_bub*Q; // scattering cross section [cm²]
  if (VERBOSE > 0) cout << "Estimated scattering cross section: " << sigma << " cm^2" << endl;
  
  // Number density of glass microspheres (very approximate)
  double V = 4.*pi/3.*(pow(r_bub,3) - pow(r_bub-d_bub,3)); // volume of glass [cm³]
  double m = rho_bub*V; // mass of hollow sphere [g]
  double N = density/m; // number concentration [cm⁻³]
  if (VERBOSE > 0) cout << "Estimated number density of bubbles: " << N << "/cm^3" << endl;
  
  // Scattering parameters
  double eps = N*sigma; // scattering coefficient [cm⁻¹]
  double L = 1./eps; // scattering length [cm]
  if (VERBOSE > 0) cout << "Estimated scattering length: " << L << " cm" << endl;

  return L;
}

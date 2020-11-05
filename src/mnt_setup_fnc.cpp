
#include <TChain.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <THStack.h>
#include <string>
#include <vector>
#include <iostream>


std::string const DATADIR = "/storage/epp2/phshgg/DVTuples__v23/";

struct path_data {
	std::string const DATADIR,
	sqrts,
	stripping,
	year,
	polarity,
	stream;
};

struct plot_config {
	std::string const expression;
	int const nbins;
	float const xmin,
	xmax;
	std::string const label,
	unit;

  //create constructor to allow struct to configure correctly in the emplace_back function
  plot_config(std::string const expression, int const nbins, float const xmin, float const xmax, std::string const label, std::string const unit)
    :expression(std::move(expression))
    ,nbins(nbins)
    ,xmin(xmin)
    ,xmax(xmax)
    ,label(std::move(label))
    ,unit(std::move(unit))
  {}
};

std::vector<plot_config> plot_configurations;


TH1F* make_TH1F(struct path_data path_in, struct plot_config hist_input, std::string hist_name, std::string hist_text) {
	TChain ch("Z/DecayTree");
	std::string path = path_in.DATADIR + path_in.sqrts + "TeV_" + path_in.year + "_" + path_in.stripping + "_" + path_in.polarity + "_" + path_in.stream + ".root";
	ch.Add(path.c_str());
	
	TH1F *hist = new TH1F(hist_name.c_str(), hist_text.c_str(), hist_input.nbins, hist_input.xmin, hist_input.xmax);
	std::string expression = hist_input.expression +">>"+ hist->GetName();
	ch.Draw(expression.c_str());
	return hist;
}

void plot_data_sim(struct path_data mnt_in, struct path_data sim_in, struct plot_config hist_input) {
	std::string hist_title = "Measurement vs Simulation for " + hist_input.label,
	  x_axis = hist_input.label + " " + hist_input.unit,
	  sim_text = hist_input.label + " Simulation Data;" + x_axis,
	  mnt_text = hist_input.label + " Experimental Data;" + x_axis;

	TH1F *hist_mnt = make_TH1F(mnt_in, hist_input, "measurement", mnt_text);
	TH1F *hist_sim = make_TH1F(sim_in, hist_input, "simulation", sim_text);
	
	hist_sim->Scale(hist_mnt->Integral()/hist_sim->Integral());
	
	TCanvas canv;
	hist_sim->SetStats(false);
	hist_sim->Draw("HIST");
	hist_mnt->Draw("SAME E");
	hist_mnt->SetLineColor(2);
	canv.BuildLegend();
	std::string const filename = "Measurement_" + hist_input.label + ".pdf";
	canv.SaveAs(filename.c_str());
}

int main() {

	/*Measurement Chain*/
	path_data const measurement_in = {DATADIR, "5", "32", "2017", "Down", "EW"};

	/*Simulation Chain*/
	path_data const simulation_in = {DATADIR, "13", "28r1", "2016", "Down", "Z_Sim09h"};

	
	//mup
	plot_configurations.emplace_back("1.e-3*mup_PT", 100, 15., 60., "mup_PT", "(GeV)");
	plot_configurations.emplace_back("mup_ETA", 100, 1.5, 5., "mup_ETA", "");
	plot_configurations.emplace_back("mup_PHI", 100, -4., 4., "mup_PHI", "");
	//mum
	plot_configurations.emplace_back("1.e-3*mum_PT", 100, 15., 60., "mum_PT", "(GeV)");
	plot_configurations.emplace_back("mum_ETA", 100, 1.5, 5., "mum_ETA", "");
	plot_configurations.emplace_back("mum_PHI", 100, -4., 4., "mum_PHI", "");
	//dimuon
	plot_configurations.emplace_back("1.e-3*Z_M", 100, 35, 120, "Z_M", "(GeV)");

	for (auto const & plot_configuration : plot_configurations) {
	  plot_data_sim(measurement_in, simulation_in, plot_configuration);
	}
	
	return 0;
}

#include "matplotlibcpp.h"
#include "pugiconfig.hpp"
#include "pugixml.hpp"
#include "pugixml.cpp"
#include <iostream>
#include <sstream>
#include <bitset>
#include <string>
#include <iomanip>
#include <string>
#include <vector>

namespace plt = matplotlibcpp;

std::vector<std::string> split(const std::string& str, int n)
{
    std::vector<std::string> substrings;
 
    for (size_t i = 0; i < str.size(); i += n) {
        substrings.push_back(str.substr(i, n));
    }
 
    return substrings;
}


double ToDecimal2kKomplement(std::string str, int bit, double ADU, double Multiplier) {
    std::string komplement; 
    std::string komplement2k; 
    std::string sign; 
    int bitCorrected = bit-1;
    double hex = 0;
    long double decimal = 0;
    sign = str.substr(0, 1) == "1" ? "-" : "+";

    if (sign == "-") {
        for(int i = 0; i < str.length(); i++) {
            komplement += str.substr(i, 1)  == "1" ? "0" : "1";
        }
        int addition = 1;
        for(int i = 0; i < komplement.length(); i++) {
            if (addition == 1) {
            if (komplement.substr(bitCorrected-i, 1)  == "1") {
                komplement2k += "0";
            } else {
                komplement2k += "1";
                addition = 0;
            }
            } else {
                komplement2k += komplement.substr(bitCorrected-i, 1);
            }
        
        }
        std::string reverse; 
            for(int i = 0; i < komplement2k.length(); i++) {
            reverse += komplement2k.substr(bitCorrected-i, 1);
        }
        komplement = reverse;
    } else {
        komplement = str;
    }

    for(int i = 0; i < bit; i++) {
        int number = std::stoi(komplement.substr(i, 1));
        int power = i+1;
        hex += number * std::pow(2, str.length()-power);
    }
 
    if (sign.compare("-") == 0) {
        decimal = -1 * hex/ADU*Multiplier;
    } else {
         decimal = hex/ADU*Multiplier;
    }
    return decimal;
}

int dataSize;

std::vector<double> TimeInterval(double Interval, int Seconds) {
    int length = Seconds/Interval;
    std::vector<double> time;

    double TimeBuffer = 0; 
    dataSize = length;
    for(int i = 0; i<length; i++){
        TimeBuffer += Interval; 
        time.push_back(TimeBuffer);
    }

    return time;
}

double CalculateTolerance(std::vector<double> set, double buffer){
    int pc = 0;
    double tol = 0;
    for (auto &point : set) {
        if (pc != 0 && pc != dataSize-1) {
            tol += std::abs(std::abs(set.at(pc))-std::abs(set.at(pc+1))); 
        }   
    pc += 1;
    }
    tol = tol/(dataSize-1)*3+buffer;

    return tol;
}

std::vector<double> SmothOut(std::vector<double> set, double tolerance) {
    int pc = 0;
    std::vector<double> smothSet;
    double tolBefore = CalculateTolerance(set, 0);
    int affected = 0;

    for (auto &point : set) {


        if (pc != 0 && pc != dataSize-1) {
            if (point < set.at(pc-1)+tolerance && point > set.at(pc-1)-tolerance && point < set.at(pc+1)+tolerance && point > set.at(pc+1)-tolerance)
            {
                smothSet.push_back((set.at(pc-1)+point+set.at(pc+1))/3);
                affected++;
            }else {
                smothSet.push_back(point);
            }
        }
        else {
                        smothSet.push_back(point);
        }
    pc += 1;
    }

    double tolAfter= CalculateTolerance(smothSet, 0);

    if (tolBefore-tolAfter > 0.00005 || tolAfter > 1) return SmothOut(smothSet, tolerance); else return smothSet;
}

void swap(double *xp, double *yp) {
    double temp = *xp;
    *xp = *yp;
    *yp = temp;
}


double BPM(std::vector<double> set){
    double value[4];
    double position[4];
    for(int i = 0; i < 4; i++ ){
    value[i] = 0;
    position[i] = 0; 
    }
    int count = 0;

    for (auto &point : set) {
    
    int smallest = 0;
    for(int i = 0; i < 4; i++){
        if (value[smallest] > value[i]) {
            smallest = i; 
        }
        }

    if (value[smallest] < point && point > set.at(count-1) && point > set.at(count+1)) {
        value[smallest] = point;
        position[smallest] = count;
    }
        count++;
    }

    int min_idx;
    for(int i = 0; i < 4-1; i++ ){
        min_idx = i;
        for(int y = i+1; y < 4; y++ ){
            if(position[y] < position[min_idx]){
                min_idx = y;
            }
            }

         if(min_idx!=1) {
            swap(&position[min_idx], &position[i]);
            swap(&value[min_idx], &value[i]);
        }
        
    }
    double time;
    for (int i = 0; i < 3; i++){
        time += position[i+1]*0.002 - position[i]*0.002;
    }
  
    return time/3*60;
    
}

int main()
{
    int bytes = 4;

    pugi::xml_document doc;

    pugi::xml_parse_result result = doc.load_file("zadatak.xml");

    double SRate = std::stod(doc.child("ECG").child("Events").child("Event").child("ECGInfo").child("SRate").attribute("Value").value());
    double ADU = std::stod(doc.child("ECG").child("Events").child("Event").child("ECGInfo").child("ADU").attribute("Value").value());
    double Multiplier = std::stod(doc.child("ECG").child("Events").child("Event").child("ECGInfo").child("Multiplier").attribute("Value").value());

    double Interval = 1/SRate; 
    std::vector<double> timeInterval = TimeInterval(Interval, 4);

    pugi::xml_node ECGData = doc.child("ECG").child("Events").child("Event").child("ECGData");

    std::vector<std::string> plots; 
    std::vector<std::string> labels; 
    std::vector<std::string> grids; 
    std::vector<std::string> xlabels; 
    std::vector<std::string> currentAxes; 


    int plot_row = 0;
    int plot_column = 0;
    double bpm = 0;

    for (pugi::xml_node_iterator it = ECGData.begin(); it != ECGData.end(); ++it)
    {
        int stop = 3;
        std::vector<std::string> tokens = split(it->attribute("Value").value() , bytes);
        int counter = 0;
        std::vector<double> data;
        std::vector<double> rawData;

        double prevData = -999;

            for (auto &token : tokens) {
                std::stringstream ss;
                ss << std::hex << token;
                unsigned n;
                ss >> n;
                std::bitset<16> b(n);                
                double currentData = ToDecimal2kKomplement(b.to_string(), 16, ADU, Multiplier);

                prevData = currentData;

                rawData.push_back(currentData) ;
                if (counter == dataSize-1) break;

                counter += 1;
            }

            std::stringstream plot;
            std::stringstream label;
            std::stringstream grid;
            std::stringstream xlabel;
            std::stringstream currentAxe;

            plot.clear();
            plot.str("");
            label.clear();
            label.str("");
            grid.clear();
            grid.str("");
            xlabel.clear();
            xlabel.str("");
            currentAxe.clear();
            currentAxe.str("");

            data = SmothOut(rawData, CalculateTolerance(rawData, 0));    
            data = SmothOut(data, CalculateTolerance(data, 0.025));
            data = SmothOut(data, CalculateTolerance(data, 0.05));
            data = SmothOut(data, CalculateTolerance(data, 0.075));
            data = SmothOut(data, CalculateTolerance(data, 0.1));

            std::string name = it->name();
            if (name == "II"){
                bpm = BPM(data);
            }

            int si = 0;
            plot << "axs[";
            plot << plot_row << ',' << plot_column;
            plot << "].plot(["; 
            for (auto it = timeInterval.begin(); it != timeInterval.end(); it++) {
                plot << timeInterval.at(si);
                if(si != dataSize) plot << ", ";
                if(si == dataSize-1) plot << "], [ ";
                si++;
            }
            si = 0;
            for (auto it = data.begin(); it != data.end(); it++) {
                if(si != dataSize) plot << data.at(si);
                if(si != dataSize) plot << ", ";
                    si++;
            }
            plot << "], c='black')";

            label << "axs[" << plot_row << ',' << plot_column << "].set_ylabel('" << it->name() << "')";
            xlabel << "axs[" << plot_row << ',' << plot_column << "].set_xlabel('t[S]')";
            grid << "axs[" << plot_row << ',' << plot_column << "].grid()";
            currentAxe << "axs[" << plot_row << ',' << plot_column << "]";




            plots.push_back(plot.str());
            labels.push_back(label.str());
            xlabels.push_back(xlabel.str());
            grids.push_back(grid.str());
            currentAxes.push_back(currentAxe.str());



            plot_row += 1; 
            if (plot_row == 6 && plot_column != 2) {plot_column += 1; plot_row = 0;}
            if(plot_row == 6 && plot_column == 2) break;

    }


Py_Initialize();
PyRun_SimpleString("import pylab");
PyRun_SimpleString("import numpy as np");


PyRun_SimpleString("px = 1/pylab.rcParams['figure.dpi']");
PyRun_SimpleString("fig, axs = pylab.subplots(figsize=(3000*px, 1400*px), nrows=6, ncols=2)");
PyRun_SimpleString("pylab.tight_layout(pad = 8, h_pad=1, w_pad=2)");
std::string python;
python = "fig.suptitle('EKG, BPM: " + std::to_string(bpm) + "',fontsize=30)";
PyRun_SimpleString(python.c_str());
PyRun_SimpleString("pylab.setp(axs, xticks=[1,2,3,4], yticks=[-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3])");

for (int i = 0; i < 12; i++) {
std::string string;
PyRun_SimpleString(plots.at(i).c_str());
PyRun_SimpleString(grids.at(i).c_str());
PyRun_SimpleString(xlabels.at(i).c_str());
PyRun_SimpleString(labels.at(i).c_str());

string = currentAxes.at(i)+".grid(color='white')";
PyRun_SimpleString(string.c_str());
string = currentAxes.at(i)+".set_facecolor('#e5e5e5')";
PyRun_SimpleString(string.c_str());
string = currentAxes.at(i)+".spines['top'].set_color('#e5e5e5')";
PyRun_SimpleString(string.c_str());
string = currentAxes.at(i)+".spines['bottom'].set_color('#e5e5e5')";
PyRun_SimpleString(string.c_str());
string = currentAxes.at(i)+".spines['left'].set_color('#e5e5e5')";
PyRun_SimpleString(string.c_str());
string = currentAxes.at(i)+".spines['right'].set_color('#e5e5e5')";
PyRun_SimpleString(string.c_str());
string = "pylab.xticks([1,2,3,4])";
PyRun_SimpleString(string.c_str());
}

PyRun_SimpleString("pylab.savefig('EKG.png', dpi=200)");

PyRun_SimpleString("pylab.show()");
Py_Exit(0);

}

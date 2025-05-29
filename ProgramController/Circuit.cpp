#include "Circuit.h"
#include <iostream>
#include <bits/stdc++.h>

using namespace std;
using namespace Controller;
using namespace Model;

double Coefficient_management(string a) {
    regex r1("(\\S*)([k|m|u|n|Meg])");
    smatch match;
    if(regex_search(a, match, r1)) {
        double output = stod(match[1].str());
        if(match[2]=="k") {
            output *=1000;
        }
        else if(match[2]=="m") {
            output *=0.001;
        }
        else if(match[2]=="u") {
            output *=0.000001;
        }
        else if(match[2]=="n") {
            output *=0.000000001;
        }
        else if(match[2]=="Meg") {
            output *=1000000;
        }
        return output;
    }
    else {
        return stod(a);
    }
}

void Circuit::addResistor(string name, string n1,string n2, string resistance) {
    try {
        Node* node1=new Node;
        node1->setName(n1);
        Node* node2=new Node;
        node2->setName(n2);
        nodes.push_back(node1);
        //nodes.push_back(node2);
        Resistor* resistor=new Resistor("R"+name,node1,node2,Coefficient_management(resistance));
        // elements.push_back(resistor);
    }catch(const exception& e) {
        cout<<e.what()<<endl;
    }
};

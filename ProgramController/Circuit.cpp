#include "Circuit.h"
#include <iostream>

using namespace std;
using namespace Controller;
using namespace Model;

void Circuit::addResistor(string name, string n1,string n2, double resistance) {
    try {
        Node* node1=new Node(n1,0);
        Node* node2=new Node(n2,0);
        Resistor* resistor=new Resistor(name,node1,node2,resistance);
        elements.push_back(resistor);
        //
    }catch(const exception& e) {
        cout<<e.what()<<endl;
    }
};

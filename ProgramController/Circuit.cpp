#include "Circuit.h"

using namespace Controller;
using namespace Model;

void Circuit::addResistor(std::string name, std::string n1, std::string n2, double resistance) {
    Node* node1 = new Node(n1, 0);
    Node* node2 = new Node(n2, 0);
    Resistor* resistor = new Resistor(name, node1, node2, resistance);
    elements.push_back(resistor);
}

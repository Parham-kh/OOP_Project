#include "Resistor.h"

using namespace Model;

Resistor::Resistor(std::string name, Node* node1, Node* node2, double resistance)
    : Element(name, "Resistor", node1, node2, 0, 0), resistance(resistance) {}

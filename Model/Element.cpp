#include "Element.h"

using namespace Model;

Element::Element(string name, string type, Node* node1, Node* node2, double voltage, double current)
    : name(name), type(type), node1(node1), node2(node2), voltage(voltage), current(current) {}

Element::~Element() = default;

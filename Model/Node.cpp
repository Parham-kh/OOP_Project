#include "Node.h"

using namespace Model;

string Node::getName() {
    return name;
}
void Node::setName(string name) {
    this->name = name;
}
double Node::getVoltage() {
    return voltage;
}
void Node::setVoltage(double voltage) {
    this->voltage = voltage;
}



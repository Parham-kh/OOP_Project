#ifndef RESISTOR_H
#define RESISTOR_H
#pragma once
#include "Element.h"
using namespace std;

namespace Model {
    class Resistor : public Element {
        double resistance = 0;

    public:
        Resistor(string name, Node* node1, Node* node2, double resistance);
    };
}

#endif //RESISTOR_H

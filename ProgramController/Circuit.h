#ifndef CIRCUIT_H
#define CIRCUIT_H
#pragma once
#include <vector>
#include <string>
#include "Model/Element.h"
#include "Model/Node.h"
#include "Model/Resistor.h"
using namespace std;
using namespace Model;

namespace Controller {
    class Circuit {
        vector<Element*> elements;
        vector<Node*> nodes;

    public:
        void addResistor(string name, string n1, string n2, double resistance);
    };
}

#endif //CIRCUIT_H

#ifndef CIRCUIT_H
#define CIRCUIT_H
#pragma once
#include <vector>
#include <string>
#include "Model/Element.h"
#include "Model/Node.h"
#include "Model/Resistor.h"

namespace Controller {
    class Circuit {
        std::vector<Model::Element*> elements;
        std::vector<Model::Node*> nodes;

    public:
        void addResistor(std::string name, std::string n1, std::string n2, double resistance);
    };
}

#endif //CIRCUIT_H

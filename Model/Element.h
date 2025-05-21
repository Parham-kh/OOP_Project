#ifndef ELEMENT_H
#define ELEMENT_H
#pragma once
#include <string>
#include "Node.h"

namespace Model {
    class Element {
        Node* node1;
        Node* node2;

    protected:
        std::string name = "";
        std::string type = "";
        double voltage = 0;
        double current = 0;

    public:
        Element(std::string name, std::string type, Node* node1, Node* node2, double voltage, double current);
        virtual ~Element();
    };
}

#endif //ELEMENT_H

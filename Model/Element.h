#ifndef ELEMENT_H
#define ELEMENT_H
#pragma once
#include <string>
#include "Node.h"
using namespace std;

namespace Model {
    class Element {
        Node* node1;
        Node* node2;

    protected:
        string name = "";
        string type = "";
        double voltage = 0;
        double current = 0;

    public:
        Element(string name, string type, Node* node1, Node* node2, double voltage, double current);
        virtual ~Element();
    };
}

#endif //ELEMENT_H

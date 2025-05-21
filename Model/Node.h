#ifndef NODE_H
#define NODE_H
#pragma once
#include <string>
using namespace std;

namespace Model {
    class Node {
        string name = "";
        double voltage = 0;

    public:
        Node(string name, double voltage);
    };
}


#endif //NODE_H

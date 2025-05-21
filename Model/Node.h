#ifndef NODE_H
#define NODE_H
#pragma once
#include <string>

namespace Model {
    class Node {
        std::string name = "";
        double voltage = 0;

    public:
        Node(std::string name, double voltage);
    };
}


#endif //NODE_H

#include <iostream>
#include "json_parser.hpp"

int main(int argc, char**argv)
{
    std::string empty;
    JsonParser parser(empty);
    parser.parse();
}

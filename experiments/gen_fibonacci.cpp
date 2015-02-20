#include <string>
#include <iostream>
#include <cstdlib>

int main(int argc, char **argv) {
    if (argc != 2) {
        std::cerr << "Usage: gen_fibonacci n." << std::endl;
        exit(1);
    }

    unsigned n = atoi(argv[1]);

    std::string s0 = "0 ";
    std::string s1 = "0 1 ";
    if (n == 0) {
        std::cout << s0 << std::endl;
    } else if (n == 1) {
        std::cout << s1 << std::endl;
    } else {
        std::string tmp;
        for (unsigned i = 2; i <= n; ++i) {
            tmp.assign(s1 + s0);
            s0.assign(s1);
            s1.assign(tmp);
        }
        std::cout << tmp << std::endl;
    }

    exit(0);
}

#include <iostream>
#include <fstream>

using namespace std;

int main()
{
	cout << "Name of file ";
	string nome;
	cin >> nome;

    ifstream abre(nome+".txt", ios::binary);
    ofstream salva(nome+".dat", ios::binary);
    float data[8];
    int numOpt;
    abre >> numOpt;
    salva.write((char*)&numOpt, sizeof(int));
    char opt;
    while(abre)
    {
        abre >> data[0] >> data[1]
             >> data[2] >> data[3]
             >> data[4] >> data[5] >> opt
             >> data[6] >> data[7];

        salva.write((char*)&data[0], sizeof(float));
        salva.write((char*)&data[1], sizeof(float));
        salva.write((char*)&data[2], sizeof(float));
        salva.write((char*)&data[3], sizeof(float));
        salva.write((char*)&data[4], sizeof(float));
        salva.write((char*)&data[5], sizeof(float));
        salva.write(&opt, sizeof(char));
        salva.write((char*)&data[6], sizeof(float));
        salva.write((char*)&data[7], sizeof(float));
    }
    return 0;
}

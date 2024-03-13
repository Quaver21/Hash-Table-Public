#include "dnadb.h"
#include <time.h>

const char BREAK[] = "*****************************************************************\n";

unsigned int hashCode(const string str){
    unsigned val = 0;
    for(unsigned i = 0; i < str.length(); i++){
        val = val * 33 + str[i];
    }
    return val;
}

string sequencer(int size){
    string sequence = "";
    for (int i = 0; i < size; i++){
        sequence += ALPHA[rand() % 4];
    }
    return sequence;
}

class Tester{
    public:
        Tester();
        int getTestCount();
        int getFailCount();
        void result(bool test);
        static bool insertTest(DnaDb dnadb, DNA arr[], int size);
        static bool removeTest(DnaDb dnadb, DNA arr[], int size);
        static bool insertRehashTest(DnaDb dnadb, DNA arr[]);
        static bool removeRehashTest(DnaDb dnadb, DNA arr[]);
        static bool getDNATest(DnaDb dnadb, DNA arr[], int size, bool answer);
        static bool hashInArray(DNA arr[], DNA dna, int size, hash_fn hash, int divisor);
        static bool inArray(DNA arr[], DNA dna, int size);
        static DnaDb copy(DnaDb& rhs);
        static DNA randDNA();
    private:
        int m_testCount;
        int m_failCount;
        static bool rehashTest(DnaDb& dnadb, DNA arr[], int size);
};

int main(){
    Tester test;
    srand(time(NULL));

    const int tableSize = 503, numDNA = tableSize / 5;
    DnaDb noCol(tableSize, hashCode), col(tableSize, hashCode), empty(tableSize, hashCode);
    DNA noColDNA[numDNA * 2], colDNA[numDNA * 2];
    
    for(int i = 0; i < numDNA / 2; i++){
        do{
            noColDNA[i] = Tester::randDNA();
        }while(Tester::hashInArray(noColDNA, noColDNA[i], i, hashCode, tableSize));
        noCol.insert(noColDNA[i]);
        do{
            colDNA[i] = Tester::randDNA();
        }while(Tester::hashInArray(colDNA, colDNA[i], i, hashCode, tableSize));
        col.insert(colDNA[i]);
    }
    for(int i = numDNA / 2; i < numDNA; i++){
        do{
            noColDNA[i] = Tester::randDNA();
        }while(Tester::hashInArray(noColDNA, noColDNA[i], i, hashCode, tableSize));
        noCol.insert(noColDNA[i]);
        do{
            colDNA[i] = DNA(colDNA[rand() % i].getSequence(), rand() % (MAXLOCID - MINLOCID + 1) + MINLOCID);
        }while(Tester::inArray(colDNA, colDNA[i], i));
        col.insert(colDNA[i]);
    }

    cout << BREAK << "Testing DnaDb::insert(DNA)\n" << BREAK << endl;
    {   cout << "Normal: Inserting non-colliding data without a rehash";
        DNA insertDNA[numDNA];
        for(int i = 0; i < numDNA; i++){
            do{
                insertDNA[i] = Tester::randDNA();
            }while(Tester::hashInArray(insertDNA, insertDNA[i], i, hashCode, tableSize)
                    || Tester::hashInArray(noColDNA, insertDNA[i], numDNA, hashCode, tableSize));
        }
        test.result(Tester::insertTest(Tester::copy(noCol), insertDNA, numDNA));
    }
    {   cout << "Edge: Inserting colliding data without a rehash";
        DNA insertDNA[numDNA];
        for(int i = 0; i < numDNA; i++){
            do{
                insertDNA[i] = DNA(colDNA[rand() % numDNA].getSequence(), rand() % (MAXLOCID - MINLOCID + 1) + MINLOCID);
            }while(Tester::inArray(insertDNA, insertDNA[i], i)
                    || Tester::inArray(colDNA, insertDNA[i], numDNA));
        }
        test.result(Tester::insertTest(Tester::copy(col), colDNA, numDNA));
    }
    {   cout << "Error: Inserting already existing data";
        test.result(Tester::insertTest(Tester::copy(noCol), noColDNA, numDNA));
    }

    cout << BREAK << "Testing DnaDb::remove(DNA)\n" << BREAK << endl;
    {   cout << "Normal: Removing non-colliding data without a rehash";
        test.result(Tester::removeTest(Tester::copy(noCol), noColDNA, numDNA * 3 / 4));
    }
    {   cout << "Edge: Removing colliding data without a rehash";
        test.result(Tester::removeTest(Tester::copy(col), colDNA, numDNA * 3 / 4));
    }
    {   cout << "Error: Removing nonexistent data";
        DNA removeDNA[numDNA];
        for(int i = 0; i < numDNA; i++){
            do{
                removeDNA[i] = Tester::randDNA();
            }while(Tester::inArray(noColDNA, removeDNA[i], numDNA));
        }
        test.result(Tester::removeTest(Tester::copy(noCol), removeDNA, numDNA));
    }

    cout << BREAK << "Testing rehashing\n" << BREAK << endl;
    {   cout << "Rehash from insertion";
        test.result(Tester::insertRehashTest(Tester::copy(noCol), noColDNA));
    }
    {   cout << "Rehash from removal";
        test.result(Tester::removeRehashTest(Tester::copy(noCol), noColDNA));
    }

    cout << BREAK << "Testing DnaDb::getDNA(string, int)\n" << BREAK << endl;
    {   cout << "Normal: Finding non-colliding data";
        test.result(Tester::getDNATest(Tester::copy(noCol), noColDNA, numDNA, true));
    }
    {   cout << "Edge: Finding colliding data";
        test.result(Tester::getDNATest(Tester::copy(col), colDNA, numDNA, true));
    }
    {   cout << "Error: Finding nonexistent data";
        DNA findDNA[numDNA];
        for(int i = 0; i < numDNA; i++){
            do{
                findDNA[i] = Tester::randDNA();
            }while(Tester::inArray(noColDNA, findDNA[i], numDNA));
        }
        test.result(Tester::getDNATest(Tester::copy(noCol), findDNA, numDNA, false));
    }

    cout << BREAK << "Number of tests: " << test.getTestCount()
         << "\nNumber of tests failed: " << test.getFailCount()
         << endl << BREAK;
}

Tester::Tester() : m_testCount(0), m_failCount(0){}

int Tester::getTestCount(){
    return m_testCount;
}

int Tester::getFailCount(){
    return m_failCount;
}

void Tester::result(bool test){
    m_testCount++;
    if(test){
        cout << "\n        TEST PASSED\n\n";
    }else{
        m_failCount++;
        cout << "\n    ****TEST FAILED\n\n";;
    }
}

bool Tester::insertTest(DnaDb dnadb, DNA arr[], int size){
    for(int i = 0; i < size; i++){
        int hash = dnadb.m_hash(arr[i].m_sequence) % dnadb.m_currentCap, correctSize = dnadb.m_currentSize, correctDeleted = dnadb.m_currNumDeleted;
        DNA* dnaptr = &dnadb.m_currentTable[hash];
        for(int j = 0; true; hash = (hash + ++j * j) % dnadb.m_currentCap){
            if(*dnaptr == EMPTY){
                if(!dnadb.insert(arr[i])
                        || !(*dnaptr == arr[i])
                        || correctSize + 1 != dnadb.m_currentSize
                        || correctDeleted != dnadb.m_currNumDeleted){
                    return false;
                }
                break;
            }else if(*dnaptr == DELETED){
                if(!dnadb.insert(arr[i])
                        || !(*dnaptr == arr[i])
                        || correctSize != dnadb.m_currentSize
                        || correctDeleted - 1 != dnadb.m_currNumDeleted){
                    return false;
                }
                break;
            }else if(*dnaptr == arr[i]){
                if(dnadb.insert(arr[i])
                        || !(*dnaptr == arr[i])
                        || correctSize != dnadb.m_currentSize
                        || correctDeleted != dnadb.m_currNumDeleted){
                    return false;
                }
                break;
            }
            dnaptr = &dnadb.m_currentTable[hash];
        }
    }
    return true;
}

bool Tester::removeTest(DnaDb dnadb, DNA arr[], int size){
    for(int i = 0; i < size; i++){
        int hash = dnadb.m_hash(arr[i].m_sequence) % dnadb.m_currentCap, correctSize = dnadb.m_currentSize, correctDeleted = dnadb.m_currNumDeleted;
        DNA* dnaptr = &dnadb.m_currentTable[hash];
        for(int j = 0; true; hash = (hash + ++j * j) % dnadb.m_currentCap){
            if(*dnaptr == EMPTY){
                if(dnadb.remove(arr[i])
                        || !(*dnaptr == EMPTY)
                        || correctSize != dnadb.m_currentSize
                        || correctDeleted != dnadb.m_currNumDeleted){
                    return false;
                }
                break;
            }else if(*dnaptr == arr[i]){
                if(!dnadb.remove(arr[i])
                        || !(*dnaptr == DELETED)
                        || correctSize != dnadb.m_currentSize
                        || correctDeleted + 1 != dnadb.m_currNumDeleted){
                    return false;
                }
                break;
            }
            dnaptr = &dnadb.m_currentTable[hash];
        }
    }
    return true;
}

bool Tester::insertRehashTest(DnaDb dnadb, DNA arr[]){
    DNA* oldAddress = dnadb.m_currentTable;
    DNA allData[dnadb.m_currentCap];
    for(int i = 0; i < dnadb.m_currentSize; i++){
        allData[i] = arr[i];
    }
    for(int i = dnadb.m_currentSize; i <= dnadb.m_currentCap / 2; i++){
        do{
            allData[i] = randDNA();
        }while(!(dnadb.getDNA(allData[i].m_sequence, allData[i].m_location) == EMPTY));
        dnadb.insert(allData[i]);
        if(dnadb.m_currentTable != oldAddress){
            return rehashTest(dnadb, allData, i);
        }
    }
    return false;
}

bool Tester::removeRehashTest(DnaDb dnadb, DNA arr[]){
    DNA* oldAddress = dnadb.m_currentTable;
    for(int i = dnadb.m_currentSize - 1; i >= dnadb.m_currentSize / 5 - 1; i--){
        dnadb.remove(arr[i]);
        if(dnadb.m_currentTable != oldAddress){
            return rehashTest(dnadb, arr, i - 1);
        }
    }
    return false;
}


bool Tester::rehashTest(DnaDb& dnadb, DNA arr[], int size){
    for(int i = 0; i < 5; i++){
        dnadb.insert(EMPTY);
    }
    if(dnadb.m_oldTable != nullptr){
        return false;
    }
    for(int i = 0; i < size; i++){
        if(dnadb.getDNA(arr[i].m_sequence, arr[i].m_location) == EMPTY){
            return false;
        }
    }
    return true;
}

bool Tester::getDNATest(DnaDb dnadb, DNA arr[], int size, bool answer){
    if(answer){
        for(int i = 0; i < size; i++){
            if(dnadb.getDNA(arr[i].m_sequence, arr[i].m_location) == EMPTY){
                return false;
            }
        }
    }else{
        for(int i = 0; i < size; i++){
            if(!(dnadb.getDNA(arr[i].m_sequence, arr[i].m_location) == EMPTY)){
                return false;
            }
        }
    }
    return true;
}

bool Tester::hashInArray(DNA arr[], DNA dna, int size, hash_fn hash, int divisor){
    for(int i = 0; i < size; i++){
        if(hash(dna.m_sequence) % divisor == hash(arr[i].m_sequence) % divisor){
            return true;
        }
    }
    return false;
}

bool Tester::inArray(DNA arr[], DNA dna, int size){
    for(int i = 0; i < size; i++){
        if(dna == arr[i]){
            return true;
        }
    }
    return false;
}

DnaDb Tester::copy(DnaDb& rhs){
    DnaDb output(rhs.m_currentCap, rhs.m_hash);
    output.m_currentSize = rhs.m_currentSize;
    output.m_currNumDeleted = rhs.m_currNumDeleted;
    output.m_oldCap = rhs.m_oldCap;
    output.m_oldSize = rhs.m_oldSize;
    output.m_oldNumDeleted = rhs.m_oldNumDeleted;
    for(int i = 0; i < output.m_currentCap; i++){
        output.m_currentTable[i] = rhs.m_currentTable[i];
    }
    if(rhs.m_oldTable != nullptr){
        output.m_oldTable = new DNA[output.m_oldCap];
        for(int i = 0; i < output.m_oldCap; i++){
            output.m_oldTable[i] = rhs.m_oldTable[i];
        }
    }
    return output;
}

DNA Tester::randDNA(){
    return DNA(sequencer(rand() % 6 + 5), rand() % (MAXLOCID - MINLOCID + 1) + MINLOCID);
}
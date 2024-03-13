/**
 * File:    dnadb.cpp
 * Project: CMSC 341 Project 4 – A DNA Database
 * Author:  Jay Buckwalter
 * Date:    5/10/22
 * Section: 03
 * E-mail:  rf29850@umbc.edu
 *
 * This file contains the implementation of the DnaDB class
 * A DnaDb is a Hash Table whose buckets each contain a DNA
 * A DNA consists of a string containing the DNA sequence and an int representing the DNA's location
 */
#include "dnadb.h"

// Name:    DnaDb (Constructor)
// Desc:    Constructor for DnaDb
// Precon:  None
// Postcon: A DnaDb will be created containing empty DNA
//          The size of the table will be the first prime number >= size within the range [MINPRIME,MAXPRIME]
DnaDb::DnaDb(int size, hash_fn hash) : m_hash(hash), m_currentSize(0), m_currentCap(findNextPrime(size - 1)),
        m_currNumDeleted(0), m_oldTable(nullptr), m_oldCap(0), m_oldSize(0), m_oldNumDeleted(0){
    m_currentTable = new DNA[m_currentCap];
}

// Name:    ~DnaDb (Destructor)
// Desc:    Destructor for DnaDb
// Precon:  None
// Postcon: All dynamically allocated memory will be deallocated
DnaDb::~DnaDb(){
    delete[] m_currentTable;
    delete[] m_oldTable;
}

// Name:    insert
// Desc:    Inserts a dna into the DnaDb
// Precon:  dna must not be within the DnaDb else doesn't insert it and returns false
// Postcon: dna will be inserted at the proper location and returns true
//          m_oldTable will be rehashed if there is an ongoing rehash or if the load factor is greater than .5
bool DnaDb::insert(DNA dna){
    bool output;
    if(dna == EMPTY || dna == DELETED){
        output = false;
    }else{ 
        output = notInOld(dna) && insertWithoutRehash(dna);
    }
    if(m_oldTable != nullptr){
        rehash25();
    }else if(lambda() > .5f){
        rehashStart();
    }
    return output;
}

// Name:    insertWithoutRehash
// Desc:    Inserts a dna into the DnaDb without rehashing
// Precon:  dna must not be within m_currentTable else does nothing and returns false
// Postcon: dna will be inserted at the proper location in m_currentTable
bool DnaDb::insertWithoutRehash(DNA& dna){
    int hash = m_hash(dna.m_sequence) % m_currentCap;
    for(int i = 0; !(m_currentTable[hash] == EMPTY || m_currentTable[hash] == DELETED); hash = (hash + ++i * i) % m_currentCap){
        if(m_currentTable[hash] == dna){
            return false;
        }
    }
    if(m_currentTable[hash] == DELETED){
        m_currNumDeleted--;
    }else{
        m_currentSize++;
    }
    m_currentTable[hash] = dna;
    return true;
}

// Name:    notInOld
// Desc:    Checks if dna is within m_oldTable
// Precon:  None
// Postcon: Returns false if dna is in m_oldTable, true otherwise
bool DnaDb::notInOld(DNA& dna){
    if(m_oldTable != nullptr){
        int hash = m_hash(dna.m_sequence) % m_oldCap;
        for(int i = 0; !(m_oldTable[hash] == EMPTY); hash = (hash + ++i * i) % m_oldCap){
            if(m_oldTable[hash] == dna){
                return false;
            }
        }
    }
    return true;
}

// Name:    remove
// Desc:    Removes a dna from the DnaDb
// Precon:  dna must be within the DnaDb else doesn't remove it and returns false
// Postcon: dna will be from the DnaDb and returns true
//          m_oldTable will be rehashed if there is an ongoing rehash
//          or if the ratio of deleted DNA to non-empty DNA is greater than .8
bool DnaDb::remove(DNA dna){
    bool output;
    if(dna == EMPTY || dna == DELETED){
        output = false;
    }else{
        output = removeFromCurrent(dna) || removeFromOld(dna);
    }
    if(m_oldTable != nullptr){
        rehash25();
    }else if(deletedRatio() > .8f){
        rehashStart();
    }
    return output;
}

// Name:    removeFromCurrent
// Desc:    Removes dna from m_currentTable
// Precon:  dna must be within m_currentTable else does nothing and returns false
// Postcon: dna will be removed from m_currentTable and returns true
bool DnaDb::removeFromCurrent(DNA& dna){
    int hash = m_hash(dna.m_sequence) % m_currentCap;
    for(int i = 0; !(m_currentTable[hash] == dna); hash = (hash + ++i * i) % m_currentCap){
        if(m_currentTable[hash] == EMPTY){
            return false;
        }
    }
    m_currNumDeleted++;
    m_currentTable[hash] = DELETED;
    return true;
}

// Name:    removeFromOld
// Desc:    Removes dna from m_oldTable
// Precon:  dna must be within m_oldTable else does nothing and returns false
// Postcon: dna will be removed from m_oldTable and returns true
bool DnaDb::removeFromOld(DNA& dna){
    if(m_oldTable == nullptr){
        return false;
    }else{
        int hash = m_hash(dna.m_sequence) % m_oldCap;
        for(int i = 0; !(m_oldTable[hash] == dna); hash = (hash + ++i * i) % m_oldCap){
            if(m_oldTable[hash] == EMPTY){
                return false;
            }
        }
        m_oldTable[hash] = DELETED;
        return true;
    }
}

// Name:    rehash25
// Desc:    Rehashes 25% of m_oldTable to m_currentTable
// Precon:  m_oldNumDeleted contains index of ongoing rehash
// Postcon: Next 25% of m_oldTable will be rehashed to m_currentTable
//          If that finishes the rehash, m_oldTable will be deleted and set to nullptr
void DnaDb::rehash25(){
    for(int i = 0; i < m_oldSize / 4 && m_oldNumDeleted < m_oldCap; m_oldNumDeleted++){
        if(!(m_oldTable[m_oldNumDeleted] == EMPTY || m_oldTable[m_oldNumDeleted] == DELETED)){
            insertWithoutRehash(m_oldTable[m_oldNumDeleted]);
            m_oldTable[m_oldNumDeleted] = DELETED;
            i++;
        }
    }
    if(m_oldNumDeleted == m_oldCap){
        delete[] m_oldTable;
        m_oldTable = nullptr;
    }
}

// Name:    rehashStart
// Desc:    Starts a rehashing of m_oldTable to m_currentTable
// Precon:  m_oldTable is nullptr
// Postcon: Sets m_oldTable to m_currentTable
//          Sets m_currentTable a new table of size newCap
//          25% of m_oldTable will be rehashed to m_currentTable
void DnaDb::rehashStart(){
    m_oldTable = m_currentTable;
    m_oldCap = m_currentCap;
    m_oldSize = m_currentSize;
    m_currentCap = findNextPrime((m_currentSize - m_currNumDeleted) * 4);
    m_oldNumDeleted = m_currentSize = m_currNumDeleted = 0;
    m_currentTable = new DNA[m_currentCap];
    rehash25();
}

// Name:    getDNA
// Desc:    Finds a DNA within the DnaDb
// Precon:  DNA must be within the DnaDb else returns EMPTY
// Postcon: Returns the DNA within DnaDb that has the m_sequence sequence and m_location location
DNA DnaDb::getDNA(string sequence, int location){
    if(m_oldTable != nullptr){
        int hash = m_hash(sequence) % m_oldCap;
        for(int i = 0; !(m_oldTable[hash] == EMPTY); hash = (hash + ++i * i) % m_oldCap){
            if(m_oldTable[hash] == DNA(sequence, location)){
                return m_oldTable[hash];
            }
        }
    }
    int hash = m_hash(sequence) % m_currentCap;
    for(int i = 0; !(m_currentTable[hash] == EMPTY || m_currentTable[hash] == DNA(sequence, location)); hash = (hash + ++i * i) % m_currentCap){}
    return m_currentTable[hash];
}

// Name:    lambda
// Desc:    Finds the load factor of m_currentTable
// Precon:  None
// Postcon: Returns the load factor of m_currentTable
float DnaDb::lambda() const{
    return (float) m_currentSize / m_currentCap;
}

// Name:    deletedRatio
// Desc:    Finds the ratio of deleted DNA in m_currentTable
// Precon:  None
// Postcon: Returns the ratio of deleted DNA to non-empty DNA in m_currentTable
float DnaDb::deletedRatio() const{
    return (float) m_currNumDeleted / m_currentSize;
}

// Name:    dump
// Desc:    Outputs a visualization of the DnaDb in array-index order
// Precon:  None
// Postcon: Contents of the DnaDb displayed to user
void DnaDb::dump() const{
    cout << "Dump for current table:\n";
    for(int i = 0; i < m_currentCap; i++){
        cout << "[" << i << "] : " << m_currentTable[i] << endl;
    }
    cout << "Dump for old table:\n";
    if(m_oldTable != nullptr)
        for(int i = 0; i < m_oldCap; i++){
            cout << "[" << i << "] : " << m_oldTable[i] << endl;
        }
}

// Name:    isPrime
// Desc:    Checks if number is prime
// Precon:  None
// Postcon: Returns true if number is prime
bool DnaDb::isPrime(int number){
    for(int i = 2; i * i <= number; i++){
        if(number % i == 0){
            return false;
        }
    }
    return true;
}

// Name:    findNextPrime
// Desc:    Finds the next prime number after i
// Precon:  i must be within [MINPRIME,MAXPRIME], else returns the closer of the two
// Postcon: Returns the smallest prime number greater than i
int DnaDb::findNextPrime(int i){
    if(i < MINPRIME){
        return MINPRIME;
    }
    else if(i >= MAXPRIME){
        return MAXPRIME;
    }else{
        while(true){
            i++;
            for(int j = 2; i % j != 0; j++){
                if(j * j > i){
                    return i;
                }
            }
        }
    }
}

// Name:    DNA (Constructor)
// Desc:    Constructor for DNA
// Precon:  None
// Postcon: A DNA will be created using the passed values
//          The location will be within [MINLOCID,MAXLOCID]
DNA::DNA(string sequence, int location){
    if((location >= MINLOCID && location <= MAXLOCID) || (location == 0 && sequence == "DELETED")){
        m_sequence = sequence;
        m_location = location;
    }else{
        m_sequence = "";
        m_location = 0;
    }
}

// Name:    getSequence
// Desc:    Gettor for m_sequence
// Precon:  None
// Postcon: Returns m_sequence
string DNA::getSequence() const{
    return m_sequence;
}

// Name:    getLocId
// Desc:    Gettor for m_location
// Precon:  None
// Postcon: Returns m_location
int DNA::getLocId() const{
    return m_location;
}

// Name:    operator= (Overloaded Assignment Operator)
// Desc:    Creates a copy of an existing DNA
// Precon:  this cannot be rhs, else nothing happens
// Postcon: this will be a deep copy of rhs
const DNA& DNA::operator=(const DNA& rhs){
    if(this != &rhs){
        m_sequence = rhs.m_sequence;
        m_location = rhs.m_location;
    }
    return *this;
}

// Name:    operator<<
// Desc:    Outputs a visualization of a DNA
//          Shows the DNA's m_sequence and m_location
// Precon:  None
// Postcon: Contents of the Post displayed to the user
ostream& operator<<(ostream& sout, const DNA &dna){
    if(!(dna == EMPTY))
        sout << dna.m_sequence << " (Location ID " << dna.m_location << ")";
  return sout;
}

// Name:    operator==
// Desc:    Compares two DNA
// Precon:  None
// Postcon: Returns true if lhs's member variables are equal to rhs's, false otherwise
bool operator==(const DNA& lhs, const DNA& rhs){
    return lhs.m_sequence == rhs.m_sequence && lhs.m_location == rhs.m_location;
}
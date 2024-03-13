CXX = g++
CXXFLAGS = -g
PROJECT = dnadb
PROJECTNAME = proj4

mytest.exe: $(PROJECT).o mytest.cpp
	$(CXX) $(CXXFLAGS) $(PROJECT).o mytest.cpp -o mytest.exe

$(PROJECT).o: $(PROJECT).h $(PROJECT).cpp
	$(CXX) $(CXXFLAGS) -c $(PROJECT).cpp

clean:
	rm *.o*
	rm *.exe
	rm *~
	rm \#*\#

run:
	./mytest.exe

val:
	valgrind ./mytest.exe

driver:
	make $(PROJECT).o
	g++ -Wall $(PROJECT).o driver.cpp -o driver.exe

test:
	valgrind ./driver.exe

submit:
	cp $(PROJECT).h $(PROJECT).cpp mytest.cpp ~/341/cs341proj/$(PROJECTNAME)

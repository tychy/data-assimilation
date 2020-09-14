kf: enzan.o kf.o
	g++ -o kf enzan.o kf.o
3henbun: enzan.o 3henbun.o
	g++ -o 3henbun enzan.o 3henbun.o
3henbun.o: 3henbun.cpp
	g++ -c 3henbun.cpp

test: enzan.o test.o
	g++ -o enzan enzan.o enzan_test.o
kf.o: kf.cpp
	g++ -c kf.cpp
enzan.o:    enzan.cpp
	g++ -c enzan.cpp
test.o:     enzan_test.cpp
	g++ -c enzan_test.cpp

clean:
	rm -f *.o enzan

echo $"Installing required softwares..."
sudo apt-get update
sudo apt-get install build-essential
sudo apt-get install cmake
sudo apt-get install libboost-dev
sudo apt-get install libglu1-mesa
sudo apt-get install libxmu-dev libxi-dev
sudo apt-get install freeglut3 freeglut3-dev
sudo apt-get install libcgal-dev

########Adding an extra boolean flag to Point_3.h
sudo sed -i '50 c\ bool flag=0;' /usr/include/CGAL/Point_3.h
cmake .
make
echo $"Now you can execute the program by running ./Peel -filename"
#include "pch.h"
#include "lab1.h"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>

//string for storing name of file to load to run program with
string filename;
//creation of vectors:
vector<JOINT> joint_vec;
vector<MEMBER> mem_vec;
vector<FORCE> ext_vec;
vector<FORCE> rxn_vec;
//creation of string for storing units used in calculations
string units;
//creation of global double pointer and single pointer Arrays to store data
double **M, **Minv, *row, *col, *E, *Answer;
int *indx, N;

using namespace std; //using namespace std to make code more readable

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Function to get Truss Size and store that information into vectors declared above.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TRUSS_SIZES getTrussSizes()
{

	string line;
	int joint_size = 0;
	int member_size = 0;
	int rxn_size = 0;
	int ext_size = 0;

	ifstream truss_size;

	while (1) {
		cout << "Name of file to load: ";
		cin >> filename;

		truss_size.open(filename);

		if (truss_size.is_open()) {
			cout << "File loaded!" << endl;
			break;
		}
		else {
			cout << "File not found, try again" << endl;
		}
	}

	getline(truss_size, line);
	line.find(JOINT_COORDINATE_HEADER);
	truss_size.ignore(256, '\n');

	while (!line.empty())
	{
		getline(truss_size, line);
		joint_size++;
	}


	getline(truss_size, line);
	line.find(MEMBER_CONNECTIVITY_HEADER);
	truss_size.ignore(256, '\n');

	while (!line.empty()) {
		getline(truss_size, line);
		member_size++;
	}

	getline(truss_size, line);
	line.find(REACTIONS_HEADER);
	truss_size.ignore(256, '\n');
	while (!line.empty()) {
		getline(truss_size, line);
		rxn_size++;
	}

	getline(truss_size, line);
	line.find(EXTERNAL_FORCES_HEADER);
	truss_size.ignore(256, '\n');
	while (!line.empty()) {
		getline(truss_size, line);
		ext_size++;
	}
	TRUSS_SIZES truss;
	truss.NJ = joint_size;
	truss.NM = member_size;
	truss.NR = rxn_size;
	truss.NEF = ext_size;

	return truss;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Function to get Truss information and data and store that information into vectors declared above.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void getTrussData(TRUSS_SIZES ts, vector<JOINT> Jf, vector<MEMBER> Mf, vector<FORCE> Ef,
	vector<FORCE> Rf, string strUnits)
{
	ifstream get_truss(filename); //create file stream from user input text file name
	string line; //string creation for pushing data into from text file
	int joint, index; //variables to allocate correct forces to joints in vectors
	int header_pos = 0; //initalize header position variable to 0
	char direction; //variable to take in reaction force direction

	//vector creation for storing headers:
	vector<string> headers;
	headers.push_back(JOINT_COORDINATE_HEADER);
	headers.push_back(MEMBER_CONNECTIVITY_HEADER);
	headers.push_back(REACTIONS_HEADER);
	headers.push_back(EXTERNAL_FORCES_HEADER);
	headers.push_back(FORCE_UNITS_HEADER);

	//While loop for going through all headers in text file
	while (header_pos < headers.size())
	{
		//stringstream temp(line);
		getline(get_truss, line);//get next line of data and store it into line
		if (line.empty())//while next line is not empty
		{
			get_truss.ignore(256, '\n');
		}
		switch (header_pos) //switch structure whith cases coresponding to header vector position
		{
		case 0:
			getline(get_truss, line);//get next line of data and store it into line
			while (!line.empty())//while next line is not empty
			{
				stringstream temp(line);
				JOINT new_joint;
				temp >> joint >> new_joint.p[0] >> new_joint.p[1];
				getline(get_truss, line);//get next line of data and store it into line
				Jf.push_back(new_joint);
			}
			joint_vec = Jf;
			header_pos++; //increment header position to find next header on next iteration of struct
			break;
		case 1:
			getline(get_truss, line);//get next line of data and store it into line
			while (!line.empty())//while next line is not empty
			{
				MEMBER new_member;
				stringstream temp1(line);
				temp1 >> joint >> new_member.j[0] >> new_member.j[1];
				getline(get_truss, line);//get next line of data and store it into line
				Mf.push_back(new_member);
			}
			mem_vec = Mf;
			header_pos++; //increment header position to find next header on next iteration of struct
			break;
		case 2:
			getline(get_truss, line);//get next line of data and store it into line
			//while loop for reaction forces vector loading
			while (!line.empty())//while next line is not empty
			{
				FORCE rxn;
				char dir;
				stringstream temp2(line);
				temp2 >> index >> rxn.j >> dir;
				if (dir == 'X')
				{
					rxn.idir = 0;
				}
				else 
				{
					rxn.idir = 1;
				}
				getline(get_truss, line);//get next line of data and store it into line
				Rf.push_back(rxn);
			}//end while loop for reaction forces vector loading
			rxn_vec = Rf; //make rxn_vec = Rf
			header_pos++; //increment header position to find next header on next iteration of struct
			break;
		case 3:
			getline(get_truss, line);//get next line of data and store it into line
			//While loop for external forces vector loading:
			while (!line.empty()) //while next line is not empty
			{
				FORCE ext; //create instane of Force struct in lab1.h
				char dir; //variable to store direction of force
				stringstream temp3(line); //creation of a stream for this sections data
				temp3 >> index >> ext.j >> ext.F >> dir; //loading of data into appropriate variables
				if (dir == 'X')
				{
					ext.idir = 0;
				}
				else {
					ext.idir = 1;
				}
				getline(get_truss, line);//get next line of data and store it into line
				Ef.push_back(ext);
			}//end while for external forces vector loading
			ext_vec = Ef;
			header_pos++; //increment header position to find next header on next iteration of struct
			break;
		case 4:
			getline(get_truss, line);//get next line of data and store it into line
			stringstream temp4(line);
			temp4 >> units; //store units string into variable
			getline(get_truss, line); //get next line of data and store it into line
			header_pos++; //increment header position to find next header on next iteration of struct
			break;
		}//end switch
	}//end while
}//end getTrussData

bool buildSystem()
{
	int rows = 0;
	float length;
	int columns = 0;
	int current_joint = 0;
	int counter = 0;
	int X_Y;
	int start_joint, end_joint;
	N = (joint_vec.size() * 2);
	E = new double[N];
	M = new double *[N];
	Minv = new double *[N];
	indx = new int[N];
	row = new double[N];
	col = new double[N];

	//creation of columns for every row in Matrix and Matrix Inverse:
	for (int i = 0; i < N; i++)
	{
		M[i] = new double[N];
		Minv[i] = new double[N];
	}

	//Building of Matrix M:
	while (rows < N)
	{

		for (X_Y = 0; X_Y <= 1; X_Y++)
		{

			for (columns; columns < mem_vec.size(); columns++) {

				if (current_joint == mem_vec[columns].j[0]) //if recorded correctly, proceed
				{
					start_joint = mem_vec[columns].j[0];
					end_joint = mem_vec[columns].j[1];
				}
				else	//if recorded reverse, flip the direction
				{
					start_joint = mem_vec[columns].j[1];
					end_joint = mem_vec[columns].j[0];
				}

				//for (columns; columns < (mem_vec.size() + rxn_vec.size()); columns++)
				//{
				if ((current_joint == start_joint) || (current_joint == end_joint))
				{
					length = sqrt((pow((joint_vec[end_joint].p[0] - joint_vec[start_joint].p[0]), 2.0) + 
						pow((joint_vec[end_joint].p[1] - joint_vec[start_joint].p[1]), 2.0)));
					M[rows][columns] = (joint_vec[end_joint].p[X_Y] - joint_vec[start_joint].p[X_Y]) / length;

				}
				else
				{
					M[rows][columns] = 0;
					//}

				}
				//cout << M[rows][columns] << endl;
			}

			for (columns; columns < mem_vec.size() + rxn_vec.size(); columns++)
			{
					if ((current_joint == rxn_vec[columns - mem_vec.size()].j) && 
						(rxn_vec[columns - mem_vec.size()].idir == X_Y)) 
					{ 
						M[rows][columns] = 1;
					}
					else {
						M[rows][columns] = 0;
					}
			}

			if (counter < ext_vec.size())
			{
				if ((X_Y == ext_vec[counter].idir) && (current_joint == ext_vec[counter].j))
				{
					E[rows] = -ext_vec[counter].F;
					counter++;
				}
				else
				{
					E[rows] = 0;
				}
			}
			else
			{
				E[rows] = 0;
			}

			columns = 0;
			rows++;
		}//end x_y for
		current_joint++; //increment current joint
	}//end building of Matrix M
	return 1;
}//end buildsystem function 

/* Function to invert matrix A with the inverse stored in Ainv.*/
void InverseNxN(double **A, int N, double **Ainv, int *indx)
{
	int i, j, k;
	double **b;

	b = new double *[N];
	for (i = 0; i < N; i++)
	{
		b[i] = new double[N];
	}


	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			b[i][j] = 0.0;
		}
	}
	for (i = 0; i < N; i++)
	{
		b[i][i] = 1.0;
	}

	ELGS(A, N, indx);

	for (i = 0; i < N - 1; i++)
	{
		for (j = i + 1; j < N; j++)
		{
			for (k = 0; k < N; k++)
			{
				b[indx[j]][k] = b[indx[j]][k] - A[indx[j]][i] * b[indx[i]][k];
			}
		}
	}

	for (i = 0; i < N; i++)
	{
		Ainv[N - 1][i] = b[indx[N - 1]][i] / A[indx[N - 1]][N - 1];
		for (j = N - 2; j >= 0; j = j - 1)
		{
			Ainv[j][i] = b[indx[j]][i];
			for (k = j + 1; k < N; k++)
			{
				Ainv[j][i] = Ainv[j][i] - A[indx[j]][k] * Ainv[k][i];
			}
			Ainv[j][i] = Ainv[j][i] / A[indx[j]][j];
		}
	}

	for (i = 0; i < N; i++)
	{
		delete[] b[i];
	}
	delete b;
}

/*******************************************************************************************************************
Function to perform the partial-pivoting Gaussian elimination. A is the original matrix in the input and transformed
matrix plus the pivoting element ratios below the diagonal in the output.  indx[] records the pivoting order.
*******************************************************************************************************************/
void ELGS(double **A, int N, int *indx)
{
	int i, j, k, itmp;
	double c1, pi, pi1, pj;
	double *c;

	c = new double[N];

	// Initialize the index
	for (i = 0; i < N; i++)
	{
		indx[i] = i;
	}

	// Find the rescaling factors, one from each row
	for (i = 0; i < N; i++)
	{
		c1 = 0;
		for (j = 0; j < N; j++)
		{
			if (fabs(A[i][j]) > c1) c1 = fabs(A[i][j]);
		}
		c[i] = c1;
	}

	// Search the pivoting (largest) element from each column 
	for (j = 0; j < N - 1; j++)
	{
		pi1 = 0;
		k = -1;
		for (i = j; i < N; i++)
		{
			pi = fabs(A[indx[i]][j]) / c[indx[i]];
			if (pi > pi1)
			{
				pi1 = pi;
				k = i;
			}
		}
		if (k == -1)
		{
			printf("Pivoting problem\n");
			exit(1);
		}

		// Interchange the rows via indx[] to record pivoting order
		itmp = indx[j];
		indx[j] = indx[k];
		indx[k] = itmp;
		for (i = j + 1; i < N; i++)
		{
			pj = A[indx[i]][j] / A[indx[j]][j];

			// Record pivoting ratios below the diagonal
			A[indx[i]][j] = pj;

			// Modify other elements accordingly
			for (k = j + 1; k < N; k++)
			{
				A[indx[i]][k] = A[indx[i]][k] - pj * A[indx[j]][k];
			}
		}
	}

	delete[] c;
}

void calculate_answer()
{
	Answer = new double [N];
	int i, j, k; //counters to step through loops of multiplying/adding
	int rows_Mat_A = N; //number of rows of matrix A
	int col_Mat_A = mem_vec.size() + rxn_vec.size();
	int rows_Mat_B = col_Mat_A; //number of columns in Matrix a -= number of rows in matrix b
	float x; //empty solution filler variable
	for (i = 0; i < rows_Mat_A; i++)
	{
		for (j = 0; j < rows_Mat_B; j++)
		{
			x = 0;
			for (k = 0; k < col_Mat_A; k++)
			{
				x = x + Minv[i][k] * E[k];
			}
			Answer[i] = x;
		}
	}
}

int main() {
	getTrussData(getTrussSizes(), joint_vec, mem_vec, ext_vec, rxn_vec, units);
	buildSystem();
	string output_name;
	cout << "What would you like to call output text file? (Exclude .txt)";
	cin >> output_name;
	output_name.append(".txt");
	ofstream outfile(output_name);
	if (!outfile)
	{
		std::cout << "oops, something went wrong! Please contact IT for support" << std::endl;
		return 0;
	}
	outfile << "Matrix M" << endl << "---------" << endl;
	printf("Matrix M:\n");
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (M[i][j] >= 0)
			{
				cout << "+";
				outfile << "+";
			}
			printf("%.2lf%c", M[i][j], (j < N - 1 ? '\t' : '\n'));//http://www.angelfire.com/oh4/psychfea/source.htm
			outfile << fixed << M[i][j] << "\t";
			
			if (j >= (N - 1))
			{
				outfile << endl;
			}
		}
	}
	InverseNxN(M, N, Minv, indx); //Creation of Inverse matrix for using to get final answer
	calculate_answer(); //Call Calculate Answer function for Matrix Multiplication 

	//Printing and Saving of Reaction and Member Forces to text file:
	printf("\n\nReaction and Member Forces:\n");
	int print_counter = 0;
	outfile << endl << endl << endl << "Solution" << endl << " ---------" << endl;
	for (int i = 0; i < N; i++)
	{
		if (i < mem_vec.size())
		{
			cout << "Member   " << i << ":  F=    ";
			outfile << "Member   " << i << ":  F=    ";
			if (Answer[i] < 0)
			{
				Answer[i] = 0 - Answer[i];
				cout << fixed << setprecision(2) << Answer[i] << " " << units << " [C]" << endl;
				outfile << fixed << setprecision(2) << Answer[i] << " " << units << " [C]" << endl;
				//printf("% .4lf%c", Answer[i], " %d", units," [C]\n");
			}
			else
			{
				cout << fixed << setprecision(2) << Answer[i] << " " << units << " [T]" << endl;
				outfile << fixed << setprecision(2) << Answer[i] << " " << units << " [T]" << endl;
				//printf("% .4lf%c", Answer[i], ' ', units,' ', '[T]', '\n');
			}
		}
		if (i >= mem_vec.size())
		{
			if (rxn_vec[print_counter].idir == 0) 
			{
				cout << "X-direction reaction on joint " << rxn_vec[print_counter].j << " =  ";
				outfile << "X-direction reaction on joint " << rxn_vec[print_counter].j << " =  ";
				if (Answer[i] >= 0)
				{
					cout << "+";
					outfile << "+";
				}
				else
				{
					//cout << "-";
					//outfile << "-";
				}
				cout << fixed << setprecision(2) << Answer[i] << " " << units << endl;
				outfile << fixed << setprecision(2) << Answer[i] << " " << units << endl;
			}
			else
			{
				cout << "Y-direction reaction on joint " << rxn_vec[print_counter].j << " =  ";
				outfile << "Y-direction reaction on joint " << rxn_vec[print_counter].j << " =  ";
				if (Answer[i] > 0)
				{
					cout << "+";
					outfile << "+";
				}
				else
				{
					cout << "-";
					outfile << "-";
				}
				cout << fixed << setprecision(2) << Answer[i] << " " << units << endl;
				outfile << fixed << setprecision(2) << Answer[i] << " " << units << endl;
			}
			print_counter++;
		}
	}
	outfile.close();
}
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <stdio.h>      /* printf */
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <math.h>       /* sqrt */
#include <iomanip>      /* std::setw */
#include <time.h>       /* time_t, time, ctime */
#include <stdio.h>      /* printf */
#include <algorithm>    // std::sort
#include <string.h>
#include <sstream>

#define INDEX(i,j,n) (i*n+j)
#define MAXINT  2147483647      /* 32-bit signed type */
#define MAXUINT 4294967295      /* 32-bit unsigned type */
#define progress
#define genseed1
//#define genseed2
#define debug
//#define cmpMatrix

using namespace std;

int vec_length = 0;
int vector_no = 0;


vector<vector<int> >
ParserPrimpoly( string poly_filename	) //read polynomial 
{
	ifstream poly_file( poly_filename.c_str(), ios::in );
	string fileline;
	if( ! poly_file ){	
		cout << "Can't open file: " << poly_filename << endl;		exit(1);
	}

	char* token;
	vector<vector<int> > primpoly;
	
	while( ! poly_file.eof() )
	{
		getline( poly_file, fileline, '\n' );

		vector<int> exponents;
		char* dup = strdup(fileline.c_str());
		token = strtok(dup, " ");
		while (token != NULL)	
		{
			exponents.push_back( atoi(token) );
			token = strtok (NULL, " ");
		}
		primpoly.push_back( exponents );
	}
	free(token);
	poly_file.close();
	return primpoly;
}

vector<string>
ParserATP( string atp_filename	) //read test cube
{
	string fileline;
	ifstream atpFile( atp_filename.c_str() , ios::in);  
	string filename;
	filename = atp_filename;

	vector<string> pattern;
	if( !atpFile )
	{  
		cout << "Can't open file: " << filename <<"\n"; 
		exit(1);	
	}
	
	getline( atpFile,fileline,'\n');
	vec_length = atoi(fileline.c_str());
	
	while( !atpFile.eof() )
	{
		getline(atpFile,fileline,'\n'); 

		if(fileline.find( "END" )==string::npos)
		{
			pattern.push_back(fileline);       // primitive vetor,  pattern<string>
		}
		else 
			break;
	}
	
	vector_no = pattern.size();
	atpFile.close();
	return pattern;

}

vector<string> DealCompModeAndModifiedPT( vector<string> pattern, vector<int> &composition_mode)
{
        vector<int> num_one, num_zero;
        /* Statistics the number of care bit of each pattern(vector)
         * to decide to integrate pattern1 and pattern2 with AND-gate(mode 1) or OR-gate(mode 0).  */
        int count1,count0;
        for (int i = 0; i < pattern.size(); i++) 
	{
                count1=0, count0=0;
                for (int j = 0; j < vec_length; j++) 
		{
                        if( pattern[i][j]=='1')                 
				count1++;
                        else if (pattern[i][j]=='0')    
				count0++;
                }
	
                if( count1 < count0 )   
			composition_mode.push_back( count1==0?0:1 );
                else    
			composition_mode.push_back( count0==0?1:0 ); //encluding count0 equal to count1.
        }
        // Statistics care-bit END

        /* Replacement the care-bit with composition_mode
         * Updated the pattern[] be a goal1[]    */
        for (int i = 0; i < pattern.size(); i++) 
	{
                for (int j = 0; j < vec_length; j++) 
		{
                        if ( composition_mode[i]==1 && pattern[i][j]=='0') 
				pattern[i][j]='2';
                        else if( composition_mode[i]==0 && pattern[i][j]=='1') 
				pattern[i][j]='2';
                }
        }

        // Replacement care-bit END
        return pattern;
}


void CalcNumCareBit( vector<string> pattern, int &cb_max, int &cb_min)
{
        cb_max=0, cb_min=MAXINT;
        for (int i = 0; i < pattern.size(); i++) 
	{
                int tmp=0;
                for (int j = 0; j < vec_length; j++) 
		{
                        if ( pattern[i][j]!='2') 
			{      
				tmp++;  
			}
                }
	
                if(cb_max<tmp)
                        cb_max = tmp;
                if(cb_min>tmp)  
			cb_min = tmp;
        }

//      printf("($INSIDE)cb_max: %2d, cb_min: %2d\n", cb_max, cb_min);

}

vector <string> Inverse( vector<string> input)
{	
	for(int i = 0 ; i < input.size() ; i++)
	{
		for(int j = 0 ; j < input[i].size() ; j++)
		{
			if(input[i][j] == '1')
				input[i][j] = '0';
			else if(input[i][j] == '0')
				input[i][j] = '1';
			else
				continue;
		}
	}

	return input;
}

bool IsPTCorrectComparator( string seed, string pattern, vector<int> xor_gate)
{
/*
        cout << "\nIn IsPTCorrectComparator():\nseed:\n";
        for (int i = 0; i < seed.length() ; i++) {
                cout << seed[i];
        }cout << endl;
        cout << "pattern:\n";
        for (int i = vec_length-1 ; i >= 0 ; i--) {
                cout << pattern[i];
        }cout << endl;
*/
//-----------------Show the data of pattern[] and seed[]. END

        string seedpattern(vec_length,'0');
        for (int i = vec_length-1; i >= 0 ; i--) {
                if ( pattern[i] != '2' && pattern[i] != seed[ seed.length()-1 ] ){
                        return false;
                }
                // LFSR ReseedingAndModifyPattern
                char c[1];
                c[0] = seed [ seed.length()-1 ];
                seedpattern[i] = seed [ seed.length()-1 ];
                for (int j = 0; j < xor_gate.size(); j++) {
                        c[0] ^= seed[ xor_gate[j] ];
                }
                for (int j = seed.length()-1; j > 0; j--) {
                        seed[j] = seed[j-1];                    }
                seed[0] = c[0];

        }
#ifdef debug
        cout << seedpattern << endl;
        cout << "----------seedpattern--------\n";
#endif

        return true;
}



void Reschedule(int *a, int *b, int m,int n)
{
        int index=0;
        int temp[m];

        int M[m][n];
        int An[m];

        for(int i=0;i<m;i++) temp[i]=i;

        for(int i=0;i<n;i++){
                int mm=index;
                int tt;
                for(int j=mm; j<m && index<m; j++){
                        if( *(a+INDEX(temp[j],i,n)) == 1 ){
                                tt=temp[index];
                                temp[index]=temp[j];
                                temp[j]=tt;
                                index++;
                        }
                }
        }

        //*************************************************
        for(int i=0;i<m;i++){
                An[i]=b[temp[i]];
                for(int j=0;j<n;j++){
                        M[i][j] = *(a+INDEX(temp[i], j, n));
                }
        }
        for(int i=0;i<m;i++){
                b[i]=An[i];
                for(int j=0;j<n;j++){
                        *(a+INDEX(i, j, n)) = M[i][j];
                }
        }
}


bool MatrixReduction(int *a,int *b, int *x,int m, int n)
{
       for (int i = 0; i < n; i++) {
                 x[i] = -1;
        }

        int pvt=0;      //pvt record the place of first 1's of this lfsr.

        for(int i=0; i<m-1; i++)
	{
                Reschedule( a+INDEX(0,0,n), b, m, n);
                //Just Reschedule the decimal value of specified_count

                for(int j=0; j<n; j++)
		{
                        if( *(a+INDEX(i, j, n) ) == 1 )
			{
                                pvt=j;
                                break;
                        }
                }

                for(int k=i+1; k<m; k++)
		{
                        //if the position [i] and [i-1] are same be 1's then...
                        if( *(a+INDEX(k,pvt,n) ) == 1 )
			{
                                for(int p=0; p<n; p++)
				{
                                        if( *(a+INDEX(i,p,n)) != *(a+INDEX(k,p,n)))
                                                *(a+INDEX(k,p,n))=1;
                                        else{ *(a+INDEX(k,p,n)) = 0 ;}
                                }
                                if(b[i]!=b[k]) 
					b[k]=1;   //b[i+1] is XORed with b[i] and b[i+1]
                                else
					b[k]=0;
                        }
                }
        }

        for(int i=m-1; i>0; i--)
	{
                pvt = -1;
                for(int j=0; j<n ;j++)
		{
                        if( *(a+INDEX(i,j,n) ) ==1 )
			{
                                pvt=j;  //pvt record the place of first 1's of this lfsr.
                                break;
                        }
                }

                if(pvt==-1) continue; //means this lfsr are all of 0's.

                for(int k=i-1 ;k>=0; k--)
		{
                        if( *(a+INDEX(k,pvt,n) ) == 1 )
			{
                                for(int p=0; p<n; p++)
				{
                                        if( *(a+INDEX(i,p,n)) != *(a+INDEX(k,p,n)))
                                                *(a+INDEX(k,p,n)) = 1;
                                        else{*(a+INDEX(k,p,n)) = 0;}
                                }
                                if(b[i]!=b[k]) b[k]=1;
                                else{b[k]=0;}
                        }
                }
        }

	int v=0;
	
        for(int i=m-1;i>=0;i--)
	{

                v=0;
                for(int j=n-1;j>=0;j--)
		{
                        if( *(a+INDEX(i,j,n)) == 1 )
                                v++; //calc the number of 1's of each seed
                }
//              if( i==4)       printf("v: %d\n", v);

                if(v==0)
		{       //no 1's bit
                        if(b[i]==0) continue;
                        else{return false;}
                }
                else if(v==1)
		{  //Just one 1's bit
                        //cout<<endl<<"1  "<<x<<"  ";
                        for(int j=n-1;j>=0;j--)
			{
                                if( *(a+INDEX(i,j,n)) == 1 )
				{

                        //              printf("x[%d]: %d, b[%d]: %d\n",j, x[j],i,b[i] );
                                        x[j]=b[i];      //record the place. of seed that bit is being 1's ,then stored
                        //              printf("x[%d]: %d, b[%d]: %d\n",j, x[j],i,b[i] );
                                        break;
                                }
                        }
                        //cout<<x<<endl;
                }
                else{   //Larger than one 1's bit
                        if(b[i]==0){
                                //cout<<endl<<"2  "<<x<<"  ";
                                for(int j=n-1;j>=0;j--){        
                                       if( *(a+INDEX(i,j,n)) ==1 )
                                                x[j]=0;  //erase same place bit in x[]
                                }
                                //cout<<x<<endl;
                        }
                        else{   //if the seed have more than one 1 bit,
                                //cout<<endl<<"3  "<<x<<"  ";
                                int pos;
                                for(int j=n-1;j>=0;j--){
                                        if(*(a+INDEX(i,j,n))==1){
                                                x[j]=0;
                                                pos=j; //then only record the leftest 1's
                                        }
                                }
                                x[pos]=1;       //then store in x[]
                                //cout<<x<<endl;
                        }
                }
        }
/*      //TEST-Method
        for (int i = 0; i < n; i++){
                if(x[i]==-1) x[i]=0; //initial val of x[] is -1, turn to 0's
                cout << x[i] << " ";
        }       cout << endl;
*/
        return true;
}

bool IsMatrixReductionSuccess(vector<vector<int> > matrix_M, string pattern, string &tseed )
{
        int a[matrix_M.size()][matrix_M[0].size()];
        int b[matrix_M.size()];

        for (int i = 0; i < matrix_M.size(); i++) 
	{
                b[i] = pattern[i]=='1'?1:0;
                for (int j = 0; j < matrix_M[i].size(); j++) 
		{
                        a[i][j] = matrix_M[i][j];
                }
        }
        int x[ matrix_M[0].size() ];

        MatrixReduction( &a[0][0] , b, x, matrix_M.size(), matrix_M[0].size() );

/*      for (int i = 0; i < x.size(); i++) {
                cout << x[i] << " ";
        }cout << endl;

*/
        for (int i = 0; i < matrix_M[0].size(); i++)
	{
                if( x[i] == -1) 
		{
                //      x[i]==1;
                        x[i] =  rand()%2 ;
                }
        }

        tseed.resize( matrix_M[0].size(), 'x' );
        for (int i = 0; i < matrix_M[0].size(); i++)
                tseed[i]= x[i]==1?'1':'0';

//      cout << "tseed: " << tseed << endl;

        return true;
}

string compressEachPTWithMatrMult( vector<vector<int> > matrix_M, string pattern  )
{
        // Initial seed value
        int minus = -1;
        vector<int> tseed( matrix_M[0].size(), minus); //seedLength

        // Erase row of matrix[] and pattern[] while mapping to  DC-set '2'
        vector<int> eraseDCStack;
        for (int i = 0; i < vec_length; i++) 
	{
                if( pattern[i] == '2' )         
			eraseDCStack.push_back( i );    
	}
	
        for (int i = 0; i < eraseDCStack.size(); i++)
	{
                matrix_M.erase( matrix_M.begin() + vec_length-1-eraseDCStack[i]  );
        }

        for (int i = vec_length-1; i >= 0;  i--) 
	{
                if( pattern[i]=='2' ) 
			pattern.erase( pattern.begin()+i );       
	}


        /* Reordering pattern[] belowing ~
         * Reordering pattern[] belowing ~
         * Reordering pattern[] belowing ~
         * */
        string tmp( pattern.length(), 'x' );
        for (int i = 0; i < pattern.length(); i++) 
	{
                tmp[ pattern.size()-i-1 ] = pattern[i];         
	}
        pattern = tmp;

#ifdef cmpMatrix
        cout << "1.Reschedule-> pattern: " << endl;
        for (int i = 0; i < pattern.length(); i++) 
	{    
		cout << pattern[i] << " ";      
	}
	cout << endl;
        printf( "matrix_M: (%3d x %3d)\n",matrix_M.size(), matrix_M[0].size() );
        MP_Reschedule_print( matrix_M , pattern);
#endif

        //20150109
        string sol;
        if(     IsMatrixReductionSuccess( matrix_M, pattern, sol) )
	{
//              cout << " Scusses !\n";
//              cout << " sol: " << sol << endl;
                return sol;
        }
	else
	{
                cout << "return no solution";
                return " $(ERROR) no solution";
        }
}



vector<string> GenSeed( vector<int> xor_gate,   int seedLength, vector<string> pattern, int &num_error, bool seed_mode  )
{

        // Create matrix_M[][]
        vector<vector<int> > matrix_M;
        matrix_M.resize( vec_length );
        
	for (int i = 0; i < vec_length; i++) 
	{
                matrix_M[i].resize( seedLength );       
	}
        
	for (int i = 0; i < seedLength; i++) 
	{
                matrix_M[i][ seedLength - 1 - i ] = 1;  
	}
	
        int h_index = 0;
        
	for (int i = seedLength; i < vec_length; i++) 
	{
                for (int j = 0; j < seedLength; j++) 
		{
                        int tmp = matrix_M[ h_index][j];
                
		        for (int k = 0; k < xor_gate.size(); k++) 
			{
                                tmp ^= matrix_M[ i-1-xor_gate[k] ] [j] ;
                        }
                        matrix_M[i][j] = tmp;
                }       
		h_index++;
        }

#ifdef genseed2
        MP_print( matrix_M , pattern[0]);
#endif
        // [note] This matrix_M can be Reuse for each original_pattern!!!!
        vector<string> seed;
        int NPass=0;
        
	for (int i = 0; i < pattern.size(); i++)        
	{
                /* Strict with size of pattern[i].size is equal to vec_length !
                 * s38417.atp always had this error.
                 */
                if( pattern[i].size() != vec_length )
		{
                        pattern[i].resize( vec_length );        
		}

                string s = pattern[i];  // Uncompress pt
                string t;                               // Compressed pt
                t  = compressEachPTWithMatrMult( matrix_M, s ) ;

#ifdef genseed1
                cout << "\n";
                cout << " ^| compressed: \n" << t << endl ;
                cout << " ^|    pattern: \n" << s  << endl ;
#endif
                seed.push_back( t );

                if ( IsPTCorrectComparator( t, s, xor_gate) )
		{
        //              printf( " Pass ! at pattern[%3d]\n", i);
                        NPass++;
                }
		else
		 {
                        if( t == " $(ERROR) no solution") 
			{
                                cout << " Error! at " << i << " $(ERROR) no solution\n";
                        }
			else if( t == "$(ERROR) no solution at MinNumColFirstly" )
			{
                        //      cout << " Error! at " << i << " with MinNumColFirstly" << endl;
                        }
			else if( t == "$(ERROR) no solution at CompressEachPT_Greedy" )
			{
                                cout    << " Error! at " << i << " with Greedy"<< endl;
                        }
                        num_error++;
                //      printf("Error at pattern %d, num_error: %d\n", i, num_error);
                        return seed;
                }
	}

        // Record the value of num_error.
        num_error = pattern.size() - NPass;

#ifdef progress
        cout << "END" <<endl;
        cout << "NumPass : " << NPass << endl;
        cout << "NumError: " << pattern.size()-NPass << endl;
        cout << "_________________________________End of GenSeed();\n\n";
#endif
	for (int i=0; i<seed.size();i++)
	{
		cout << seed[i] << endl;
	}
        return seed;
}

vector<string> ReseedingAndModifyPattern( vector<string> pattern, vector<string> original_pattern, vector<string> seed, int lfsr_size, vector<int> xor_gate  )
{
        vector<string> seed_pattern;

        for(int i=0;i<vector_no;i++){
                string temp_seed(vec_length,'0');

                int lfsr_size = seed[i].length();
                bool flag=false;
                string zero_seed(vec_length,'0');       //initial
                for(int j = vec_length-1; j>=0; j--){
                        char c[1];
                        temp_seed[j] = seed[i][ lfsr_size-1 ];
                        c[0] = seed[i][ lfsr_size-1 ];
                        for (int k = 0; k < xor_gate.size(); k++) {
                                c[0] ^= seed[i][ xor_gate[k] ];
                        }
                        for(int l = lfsr_size-1; l>0; l--)
                                seed[i][l]=seed[i][l-1];
                        seed[i][0] = c[0];
                }
                seed_pattern.push_back(temp_seed);
        }
	
     return seed_pattern;
}

vector< string > CalFinalResult(vector<string> input1, vector<string> input2, vector<int> compositionMode)
{
	vector<string>result;
	

	for(int i = 0 ; i < input1.size() ; i++)
	{
		string temp = "";	
		for(int j = 0 ; j < input1[i].size() ; j++)
		{
			stringstream ss1;
			stringstream ss2;
			stringstream ss3;
			int tmp1 ;
			int tmp2 ;	
			int tmp3 ;
			ss1 << input1[i][j];
			ss1 >> tmp1;
			ss2 << input2[i][j];
			ss2 >> tmp2;
			if(compositionMode[i] == 0)
			{
				tmp3 = tmp1 | tmp2;
			}			
			else
			{
				tmp3 = tmp1 & tmp2;				
			}
			char tmp4 ;
			ss3 << tmp3;
			ss3 >> tmp4;
			temp += tmp4;
		}
		result.push_back(temp);
	}

/*	
	for(int i = 0 ; i < result.size() ; i++)
		cout<<result[i]<<endl;
*/	
	return result;	
}

vector<int> Compare( vector<string> originalPattern, vector<string>result1 , vector<string>result2,  vector<string>result3,  vector<string>result4,int cbMax)
{
	vector<int>matchResult;
	ofstream fout("transtion1",ios::out);
	ofstream fout2("length" ,ios::out);
	fout2 << originalPattern[1].size();
	
	for(int i = 0 ; i < originalPattern.size() ; i++)
	{
		bool success1 = true;
		bool success2 = true;
		bool success3 = true;
		bool success4 = true;
		int count1 = 0;
 		int count2 = 0;
		int count3 = 0;
		int count4 = 0;
		int count5 = 0;
		for(int j = 0 ; j < originalPattern[i].size() ; j++)
		{
			if(originalPattern[i][j] == '2')
				continue;
			else if ( originalPattern[i][j] == '1')
			{
				if(result1[i][j] == '1')
					continue;
				else 
				{
					success1 = false;
					break;
				}
			}
			else
			{
				if(result1[i][j] == '0')
					continue;
				else 
				{
					success1 = false;
					break;
				}
			}
			
			if(success1 == false)
				break;
		}
	
		for(int j = 0 ; j < originalPattern[i].size() ; j++)
		{
			if(originalPattern[i][j] == '2')
				continue;
			else if ( originalPattern[i][j] == '1')
			{
				if(result2[i][j] == '1')
					continue;
				else 
				{
					success2 = false;
					break;
				}
			}
			else
			{
				if(result2[i][j] == '0')
					continue;
				else 
				{
					success2 = false;
					break;
				}
			}
			
			if(success2 == false)
				break;
		}

		for(int j = 0 ; j < originalPattern[i].size() ; j++)
		{
			if(originalPattern[i][j] == '2')
				continue;
			else if ( originalPattern[i][j] == '1')
			{
				if(result3[i][j] == '1')
					continue;
				else 
				{
					success3 = false;
					break;
				}
			}
			else
			{
				if(result3[i][j] == '0')
					continue;
				else 
				{
					success3 = false;
					break;
				}
			}
			
			if(success3 == false)
				break;
		}
	
		for(int j = 0 ; j < originalPattern[i].size() ; j++)
		{
			if(originalPattern[i][j] == '2')
				continue;
			else if ( originalPattern[i][j] == '1')
			{
				if(result4[i][j] == '1')
					continue;
				else 
				{
					success4 = false;
					break;
				}
			}
			else
			{
				if(result4[i][j] == '0')
					continue;
				else 
				{
					success4 = false;
					break;
				}
			}
			
			if(success4 == false)
				break;
		}
		if( (success1 || success2 || success3 || success4) == true)
		{
			int one = 1;
			matchResult.push_back(one);
		}
		else
		{
			int zero = 0;
			matchResult.push_back(zero);
		}
		if(success1 == true)
		{
			for(int j = 0 ; j < originalPattern[i].size()-1; j++)
			{
				if(result1[i][j] != result1[i][j+1])
					count1++;
			}
		}
		if(success2 == true)
		{
			for(int j = 0 ; j < originalPattern[i].size()-1; j++)
			{
				if(result2[i][j] != result2[i][j+1])
					count2++;
			}
		}

		if(success3 == true)
		{
			for(int j = 0 ; j < originalPattern[i].size()-1; j++)
			{
				if(result3[i][j] != result3[i][j+1])
					count3++;
			}
		}
		if(success4 == true)
		{
			for(int j = 0 ; j < originalPattern[i].size()-1; j++)
			{
				if(result4[i][j] != result4[i][j+1])
					count4++;
			}
			/*
			cout << "result4:" ;
			for(int j = 0 ; j < originalPattern[i].size()-1; j++)
			{
				cout << result4[i][j];
			}*/
		}

		int best = MAXINT;
		if(count1 != 0 && best > count1)
			best = count1;
		if(count2 != 0 && best > count2)
			best = count2;
		if(count3 != 0 && best > count3)
			best = count3;
		if(count4 != 0 && best > count4)
			best = count4;
		
		int out = 0;
		if(best == count1 )
			out = 1;
		else if (best == count2)
			out = 2;
		else if (best == count3)
			out = 3;
		else if (best == count4)
			out = 4;
		if(out == 1)
		{	
			/*for(int j = 0 ; j < originalPattern[i].size() ; j++)
			{
				fout << result1[i][j];
			}*/
			fout << cbMax << ' ';
			fout << count1;
			fout << endl;
		}
		else if(out == 2)
		{
			/*for(int j = 0 ; j < originalPattern[i].size() ; j++)
			{
				fout << result2[i][j];
			}*/
			fout << cbMax << ' ';
			fout << count2;
			fout << endl;
		}
		else if(out == 3)
		{	
			/*for(int j = 0 ; j < originalPattern[i].size() ; j++)
			{
				fout << result3[i][j];
			}*/
			fout << cbMax << ' ';
			fout << count3;
			fout << endl;
		}
		else if(out == 4)
		{
			/*for(int j = 0 ; j < originalPattern[i].size() ; j++)
			{
				fout << result4[i][j];
			}*/
			fout << cbMax << ' ';
			fout << count4;
			fout << endl;
		}
		else
		{
			for(int j = 0 ; j < originalPattern[i].size() ; j++)
			{
				if(originalPattern[i][0]=='2')
					originalPattern[i][j]='0';
				if(originalPattern[i][j]=='2')
					originalPattern[i][j]=originalPattern[i][j-1];
			//	fout << originalPattern[i][j];
			}
			for(int j = 0 ; j < originalPattern[i].size()-1; j++)
                        {
                                if(originalPattern[i][j] != originalPattern[i][j+1])
                                        count5++;
                        }
			fout << "-1 ";
			fout << count5;
			fout << endl;
		}
	}/*
	for(int j = 0 ; j < matchResult.size(); j ++)
	{
		cout << matchResult[j];
	}
	cout<<endl;*/
	fout.close();
	return matchResult;
}		

int main( int argc, char *argv[])
{	
	if( argc <= 1)
	{
		cout << "Error! Invalid parameter!" << endl;
		return 0;
	}

	string atpFileName = argv[1];
	string polyFileName = "min-term-168.poly";
	string polyFileName2 = "three-term-1656.poly";
	string polyFileName3 = "five-term-1656.poly";

	clock_t start_tick, end_tick; 
	double elapsed; 
	start_tick = clock();  	
	
	vector< vector<int> > primPoly, primPoly2, primPoly3;	
	vector< string > originalPattern, maskedPattern, maskedInversePattern;
	vector< string > reSeedingPattern1, reSeedingPattern2, reSeedingInversePattern1, reSeedingInversePattern2;
	vector< string > result1, result2, result3, result4;
	vector< int > compositionMode;
	vector< int > xorGate1, xorGate2;
	vector< string > seed1, seed2, inverseSeed1, inverseSeed2;
	int cbMax , cbMin ;

	originalPattern = ParserATP( atpFileName );
	primPoly = ParserPrimpoly( polyFileName );
	primPoly2 = ParserPrimpoly( polyFileName2 );
	primPoly3 = ParserPrimpoly( polyFileName3 );
/*
	cout << "vec_length is " << vec_length << endl;
	for(int i = 0 ; i < originalPattern.size() ; i++)
		cout << " pattern is " << originalPattern[i] << endl;	
*/	
/*
	for(int i = 0 ; i < primPoly.size(); i++)
	{
		for(int j = 0 ; j < primPoly[i].size(); j++)
			cout << primPoly[i][j] << " " ;
		cout << endl;
	}
*/
	maskedPattern = DealCompModeAndModifiedPT( originalPattern, compositionMode );
	maskedInversePattern = Inverse( maskedPattern );
/*	
	for(int i = 0 ; i < maskedPattern.size() ; i++)
		cout << " maskedpattern is " << maskedPattern[i] << endl;	
*/
/*	for(int i = 0 ; i < maskedInversePattern.size(); i ++)
		cout << "maskedInversePattern is " << maskedInversePattern[i] << endl;
*/
	CalcNumCareBit(maskedPattern, cbMax, cbMin );

	int numError;
	int lengLFSR = 0;
	bool done = false;


	
	for(int i = cbMax ; i < primPoly3.size(); i++)
	{
		numError = 0;
		
		int length1 = primPoly3[i-1][0];
		xorGate1 = primPoly3[i-1];
		xorGate1.erase( xorGate1.begin() ); 
		
		if( length1 > vec_length )
			break;

		seed1 = GenSeed( xorGate1, length1, maskedPattern, numError, 0);	
/*
		for(int i = 0 ; i < seed1.size(); i++)
			cout << "seed1 : " << seed1[i] << endl;
*/
		if(numError > 0)
			continue;		
		
		reSeedingPattern1 = ReseedingAndModifyPattern( maskedPattern, originalPattern, seed1, length1, xorGate1);
		
		if(numError == 0)
			break;
	}
/*
	for(int i = 0 ; i < reSeedingPattern1.size() ; i++)
		cout<<"reSeedingPattern1["<<i<<"]"<<reSeedingPattern1[i]<<endl;
	cout<<endl;
*/
	for(int i = cbMax ; i < primPoly3.size(); i++)
	{
		numError = 0;
		
		int length1 = primPoly3[i-1][0];
		xorGate1 = primPoly3[i-1];
		xorGate1.erase( xorGate1.begin() ); 
		
		if( length1 > vec_length )
			break;

		inverseSeed1 = GenSeed( xorGate1, length1, maskedInversePattern, numError, 0);	
/*
		for(int i = 0 ; i < inverseSeed1.size(); i++)
			cout << "inverseSeed1 : " << inverseSeed1[i] << endl;
*/
		if(numError > 0)
			continue;
		
		reSeedingInversePattern1 = ReseedingAndModifyPattern(maskedInversePattern, originalPattern, inverseSeed1, length1, xorGate1);

		if(numError == 0)
			break;
	}
	
	reSeedingInversePattern1 = Inverse(reSeedingInversePattern1);
/*
	for(int i = 0 ; i < reSeedingInversePattern1.size() ; i++)
		cout<<"reSeedingInversePattern1["<<i<<"]"<<reSeedingInversePattern1[i]<<endl;
	cout<<endl;
*/
	for(int i = cbMax ; i < primPoly2.size(); i++)
	{
		numError = 0;
		
		int length1 = primPoly2[i-1][0];
		xorGate1 = primPoly2[i-1];
		xorGate1.erase( xorGate1.begin() ); 
		
		if( length1 > vec_length )
			break;

		seed2 = GenSeed( xorGate1, length1, maskedPattern, numError, 0);	
/*
		for(int i = 0 ; i < seed2.size(); i++)
			cout << "seed2 : " << seed2[i] << endl;
*/
		if(numError > 0)
			continue;		

		reSeedingPattern2 = ReseedingAndModifyPattern(maskedPattern, originalPattern, seed2, length1, xorGate1);

		if(numError == 0)
			break;
	}
/*
	for(int i = 0 ; i < reSeedingPattern2.size() ; i++)
		cout<<"reSeedingPattern2["<<i<<"]"<<reSeedingPattern2[i]<<endl;
	cout<<endl;
*/	
	for(int i = cbMax ; i < primPoly2.size(); i++)
	{
		numError = 0;
		
		int length1 = primPoly2[i-1][0];
		xorGate1 = primPoly2[i-1];
		xorGate1.erase( xorGate1.begin() ); 
		
		if( length1 > vec_length )
			break;

		inverseSeed2 = GenSeed( xorGate1, length1, maskedInversePattern, numError, 0);	
/*
		for(int i = 0 ; i < inverseSeed2.size(); i++)
			cout << "inverseSeed2 : " << inverseSeed2[i] << endl;
*/
		if(numError > 0)
			continue;		

		reSeedingInversePattern2 = ReseedingAndModifyPattern(maskedInversePattern, originalPattern, inverseSeed2, length1, xorGate1);

		if(numError == 0)
			break;
	}
		reSeedingInversePattern2 = Inverse(reSeedingInversePattern2);
/*
	for(int i = 0 ; i < reSeedingInversePattern2.size() ; i++)
		cout<<"reSeedingInversePattern2["<<i<<"]"<<reSeedingInversePattern2[i]<<endl;
	cout<<endl;
*/	

	result1 = CalFinalResult(reSeedingPattern1, reSeedingPattern2, compositionMode);
	result2 = CalFinalResult(reSeedingInversePattern1, reSeedingPattern2, compositionMode);
	result3 = CalFinalResult(reSeedingInversePattern1, reSeedingInversePattern2, compositionMode);
	result4 = CalFinalResult(reSeedingPattern1, reSeedingInversePattern2, compositionMode);

	vector <int> compareResult;
	compareResult =  Compare(originalPattern , result1 , result2 , result3 , result4,cbMax);

/*	
	for(int i = 0 ; i < compareResult.size(); i ++)
	{
		cout << i << " : " <<compareResult[i]<<endl;
	}
*/
	
			
}

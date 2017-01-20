package cppcode;

use Inline CPP => config => namespace => 'cppcode';
use Inline CPP => <<'END_CPP';

using namespace std;
#include <string>

/* last checked Jan2017 */
/* B Contreras-Moreira EEAD-CSIC */

class Consensus
{
	public:
		Consensus() : 
			name(""),
			source1(""),source2(""),
			str1(""),str2(""),
			consensus(""),evidence(""),
			min_overlap(1),
			priority(2) {}
			
		~Consensus(){}

  		void set_priority(unsigned int prior)
    	{
        	priority = prior;
    	}
	
		void set_min_overlap(unsigned int overlap)
    	{
        	min_overlap = overlap;
    	}
			
		void set_sources(char *label1, char *label2)
    	{
        	source1 = label1;
			  source2 = label2;
    	}
	
		void set_sequences(char *seqname, char *seq1, char *seq2)
    	{
        	name = seqname; 
			str1 = seq1;
			str2 = seq2;
    	}
		
		const char *get_consensus()
    	{
        	return consensus.c_str();
    	}
		
		const char *get_evidence()
    	{
        	return evidence.c_str();
    	}
		
		int calc_consensus()
  		{
			int len1 = str1.size(); 
			int len2 = str2.size(); 
			int length=0, align[5];
			int *curr = new int [len2];
			int *prev = new int [len2];
			int *swap = NULL;
      
			
      if(len1==0 || len2==0) return 0;
		     
			for(int i = 0; i<len1; ++i)
			{
				for(int j = 0; j<len2; ++j)
				{
					if(str1[i] == str2[j] || str2[j] == 'N')
					{
						if(i == 0 || j == 0) curr[j] = 1;
						else curr[j] = 1 + prev[j-1];
		                               
						if(length < curr[j])
						{
							length = curr[j];
							align[0] = length;
							align[1] = i-length+1;
							align[2] = i;
							align[3] = j-length+1;
							align[4] = j;
						}
					}
					else curr[j] = 0;
				}
				swap=curr;
				curr=prev;
				prev=swap;
			}	
			     
			delete [] curr;
			delete [] prev;
       
			// construct consensus based on alignment
			--len1; // convert to 0-based
			--len2;
			if(align[0] >= min_overlap)
			{
				if(align[2] == len1 && align[3] == 0)
				{
					// 1----------
					//      2-----------
					consensus = str1.substr(0,align[1]) + str2;	
					evidence  = source1 + '.' + source2;	
				}
				else if(align[1] == 0 && align[4] == len2)
				{
					//     1----------
					// 2-----------
					consensus = str2 + str1.substr(align[2]+1,len1-align[2]);
					evidence  = source2 + '.' + source1;
				}
				else if(align[3] == 0 && align[4] == len2)
				{
					// 1--------------------
					//      2-----------
					consensus = str1;
					evidence  = source1 + '<' + source2;		
				}
				else if(align[1] == 0 && align[2] == len1)
				{
					//     1----------
					// 2------------------
					consensus = str2;
					evidence  = source2 + '<' + source1;
				}
				else
				{
					if(priority == 1)
					{ 
						consensus = str1;
						evidence  = source1 + "-mismatches";
					}	
					else 
					{
						consensus = str2;
						evidence  = source2 + "-mismatches";
					}	
				} 
			}
			else
			{
				if(priority == 1)
				{ 
					consensus = str1;
					evidence  = source1 + "-noover";
				}	
				else 
				{
					consensus = str2;
					evidence  = source2 + "-noover";
				}
			}     
			
			return length;
    	}

	private:
		string name; 
		string source1;
		string source2;
		string str1;
		string str2;
		string consensus;
		string evidence;
		unsigned int min_overlap;
		unsigned int priority;	
};



END_CPP

1;

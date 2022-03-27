using System;
using System.Collections.Generic;
using System.Text;
using System.Linq;

namespace Lab2
{
    public class DnaMatching : MarshalByRefObject
    {
        /// <summary>
        ///   Wariant I z prostym systemem oceny jakości dopasowania dwóch sekwencji DNA
        /// </summary>
        /// <param name="seq1"> pierwsza niepusta sekwencja DNA złożona ze znaków 'A', 'C', 'G', 'T'</param>
        /// <param name="seq2"> druga niepusta sekwencja DNA złożona ze znaków 'A', 'C', 'G', 'T'</param>
        /// <returns>(dopasowanie [ciąg 1], dopasowanie [ciąg 2], wartość całego dopasowania). 
        ///  w pierwszym etapie można zwracać nulle zamiast ciągów dopasowania </returns>
        public (string matchingSeq1, string matchingSeq2, int bestMatchingValue) FindMatchingV1(string seq1, string seq2)
        {
            const int matchValue = 1;
            const int mismatchValue = -3;
            const int gapValue = -2;

            // Metoda do zaimplementowania w etapach 1 i 2
            seq1 = " " + seq1;
            seq2 = " " + seq2;

            int[,] tab = new int[seq1.Length, seq2.Length];
            string[,] pathTab = new string[seq1.Length, seq2.Length];

            int i;
            int j;
            for(i = 0; i < seq1.Length; i++)
            {
                for(j = 0; j < seq2.Length; j++)
                {
                    if (i == 0 && j == 0) { tab[i, j] = 0; pathTab[i, j] = "start"; continue; }// wypelnienie lewego gornego rogu (0,0)
                    else if (i == 0) { tab[i, j] = tab[i, j - 1] + gapValue; pathTab[i, j] = "l"; continue; }// wypelnianie pierwszego wiersza
                    else if (j == 0) { tab[i, j] = tab[i - 1, j] + gapValue; pathTab[i, j] = "g"; continue; }// wypelnianie pierwszej kolumny


                    int bestVal;
                    string bestMove;

                    // idziemy w prawo
                    bestVal = tab[i, j -1] + gapValue;
                    bestMove = "l";
                    // idziemy w prawo-dol
                    if(seq1[i] == seq2[j])
                    {
                        if (tab[i - 1, j - 1] + matchValue > bestVal)
                        {
                            bestVal = tab[i - 1, j - 1] + matchValue;
                            bestMove = "lg";
                        }
                    }
                    else
                    {
                        if (tab[i - 1, j - 1] + mismatchValue > bestVal)
                        {
                            bestVal = tab[i - 1, j - 1] + mismatchValue;
                            bestMove = "lg";
                        }
                    }
                    // idziemy w dol
                    if(tab[i - 1, j] + gapValue > bestVal)
                    {
                        bestVal = tab[i - 1, j] + gapValue;
                        bestMove = "g";
                    }

                    tab[i, j] = bestVal;
                    pathTab[i, j] = bestMove;
                }
            }

            // odtwarzanie sekwencji ze sciezki
            i = seq1.Length - 1;
            j = seq2.Length - 1;
            string res_seq1 = "";
            string res_seq2 = "";
            while(pathTab[i, j] != "start")
            {
                if(pathTab[i, j] == "g") //idziemy do gory
                {
                    res_seq1 = seq1[i] + res_seq1;
                    res_seq2 = "-" + res_seq2;
                    i--;
                }
                else if(pathTab[i, j] == "lg") // idziemy lewo gora
                {
                    res_seq1 = seq1[i] + res_seq1;
                    res_seq2 = seq2[j] + res_seq2;
                    i--;
                    j--;
                }
                else if(pathTab[i, j] == "l") // dzieiym w lewo
                {
                    res_seq1 = "-" + res_seq1;
                    res_seq2 = seq2[j] + res_seq2;
                    j--;
                }
            }
            return (res_seq1, res_seq2, tab[seq1.Length - 1, seq2.Length - 1]);
        }


        /// <summary>
        ///   Wariant II z zaawansowanym systemem oceny jakości dopasowania dwóch sekwencji DNA
        /// </summary>
        /// <param name="seq1"> pierwsza niepusta sekwencja DNA złożona ze znaków 'A', 'C', 'G', 'T'</param>
        /// <param name="seq2"> druga niepusta sekwencja DNA złożona ze znaków 'A', 'C', 'G', 'T'</param>
        /// <returns>(dopasowanie [ciąg 1], dopasowanie [ciąg 2], wartość całego dopasowania). 
        ///  w trzecim etapie można zwracać nulle zamiast ciągów dopasowania </returns>
        public (string matchingSeq1, string matchingSeq2, int bestMatchingValue) FindMatchingV2(string seq1, string seq2)
        {
            const int matchValue = 1;
            const int mismatchValue = -3;
            const int gapStartValue = -5;
            const int gapContinuationValue = -2;

            // Metoda do zaimplementowania w etapach 3 i 4
            seq1 = " " + seq1;
            seq2 = " " + seq2;

            int[,] tabDown = new int[seq1.Length, seq2.Length];            
            int[,] tabRight = new int[seq1.Length, seq2.Length];
            int[,] tab = new int[seq1.Length, seq2.Length];
            
            string[,] pathDownTab = new string[seq1.Length, seq2.Length];
            string[,] pathRightTab = new string[seq1.Length, seq2.Length];
            string[,] pathTab = new string[seq1.Length, seq2.Length];

            int i, j;

            tab[0, 0] = 0;
            tabDown[0, 0] = -2_000_000;
            tabRight[0, 0] = -2_000_000;

            pathDownTab[0, 0] = "start";
            pathRightTab[0, 0] = "start";
            pathTab[0, 0] = "start";

            for (i = 1; i < seq1.Length; i++)
            {
                tab[i, 0] = gapStartValue + (i - 1) * gapContinuationValue;
                tabDown[i, 0] = gapStartValue + (i - 1) * gapContinuationValue;
                tabRight[i, 0] = -2_000_000;

                pathDownTab[i, 0] = "kolumna";
                pathRightTab[i, 0] = "kolumna";
                pathTab[i, 0] = "kolumna";
            }

            for(j = 1; j < seq2.Length; j++)
            {
                tab[0, j] = gapStartValue + (j - 1) * gapContinuationValue;
                tabDown[0, j] = -2_000_000;
                tabRight[0, j] = gapStartValue + (j - 1) * gapContinuationValue;

                pathDownTab[0, j] = "wiersz";
                pathRightTab[0, j] = "wiersz";
                pathTab[0, j] = "wiersz";
            }


            for (i = 1; i < seq1.Length; i++)
            {
                for (j = 1; j < seq2.Length; j++)
                {                 
                    //tabRight[i, j] = Math.Max( tab[i, j - 1] + gapStartValue, tabRight[i, j - 1] + gapContinuationValue );
                    if(tab[i, j - 1] + gapStartValue > tabRight[i, j - 1] + gapContinuationValue)
                    {
                        tabRight[i, j] = tab[i, j - 1] + gapStartValue;
                        pathRightTab[i, j] = "tab";    
                    }
                    else
                    {
                        tabRight[i, j] = tabRight[i, j - 1] + gapContinuationValue;
                        pathRightTab[i, j] = "tabRight";
                    }

                    //tabDown[i, j]  = Math.Max( tab[i - 1, j] + gapStartValue, tabDown[i - 1, j] + gapContinuationValue );
                    if(tab[i - 1, j] + gapStartValue > tabDown[i - 1, j] + gapContinuationValue)
                    {
                        tabDown[i, j] = tab[i - 1, j] + gapStartValue;
                        pathDownTab[i, j] = "tab";
                    }
                    else
                    {
                        tabDown[i, j] = tabDown[i - 1, j] + gapContinuationValue;
                        pathDownTab[i, j] = "tabDown";
                    }
                    int currMatchVal;
                    if (seq1[i] == seq2[j]) currMatchVal = matchValue;
                    else currMatchVal = mismatchValue;

                    //tab[i, j] = Math.Max(tab[i - 1, j - 1] + currMatchVal, Math.Max(tabDown[i - 1, j - 1] + currMatchVal, tabRight[i - 1, j - 1] + currMatchVal));

                    // prawo-dol
                    int bestVal = tab[i - 1, j - 1] + currMatchVal;
                    pathTab[i, j] = "tab";

                    if (tabDown[i - 1, j - 1] + currMatchVal > bestVal)
                    {
                        bestVal = tabDown[i - 1, j - 1] + currMatchVal;
                        pathTab[i, j] = "tabDown";
                    }

                    if (tabRight[i - 1, j - 1] + currMatchVal > bestVal)
                    {
                        bestVal = tabRight[i - 1, j - 1] + currMatchVal;
                        pathTab[i, j] = "tabRight";
                    }

                    tab[i, j] = bestVal;
                }
            }
            string finish;
            int result;
            result = tab[seq1.Length - 1, seq2.Length - 1];
            finish = "tab";

            if(tabDown[seq1.Length - 1, seq2.Length - 1] > result)
            {
                result = tabDown[seq1.Length - 1, seq2.Length - 1];
                finish = "tabDown";
            }
            if(tabRight[seq1.Length - 1, seq2.Length - 1] > result)
            {
                result = tabRight[seq1.Length - 1, seq2.Length - 1];
                finish = "tabRight";
            }    

             

            // odtwarzanie sekwencji ze sciezki
            i = seq1.Length - 1;
            j = seq2.Length - 1;
            string res_seq1 = "";
            string res_seq2 = "";
            string currentSet = finish;
            while ((i > 0 || j > 0) && currentSet != "start")
            {
                
                if(currentSet == "tabRight")
                {
                    res_seq1 = "-" + res_seq1;
                    res_seq2 = seq2[j] + res_seq2;

                    currentSet = pathRightTab[i, j--];
                }
                else if(currentSet == "tabDown")
                {
                    res_seq1 = seq1[i] + res_seq1;
                    res_seq2 = "-" + res_seq2;

                    currentSet = pathDownTab[i--, j];
                }
                else if(currentSet == "tab")
                {
                    res_seq1 = seq1[i] + res_seq1;
                    res_seq2 = seq2[j] + res_seq2;

                    currentSet = pathTab[i--, j--];                    
                }
                else if(currentSet == "kolumna")
                {
                    res_seq1 = seq1[i] + res_seq1;
                    res_seq2 = "-" + res_seq2;

                    currentSet = pathTab[i--, 0];
                }

                else if(currentSet == "wiersz")
                {
                    res_seq1 = "-" + res_seq1;
                    res_seq2 = seq2[j] + res_seq2;

                    currentSet = pathTab[0, j--];
                }

            }
            res_seq1.Trim();
            res_seq2.Trim();
            seq1 = "";
            seq2 = "";
            
            for(i = 0; i < res_seq1.Length; i++)
            {
                if(res_seq1[i] ==' ')
                {
                    seq1 += "-";
                }
                else
                {
                    seq1 += res_seq1[i];
                }
            }
            for (j = 0; j < res_seq2.Length; j++)
            {
                if (res_seq2[j] == ' ')
                {
                    seq2 += "-";
                }
                else
                {
                    seq2 += res_seq2[j];
                }
            }


            return (seq1, seq2, result);
        }
    }
}
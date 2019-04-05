{\rtf1\ansi\ansicpg1252\cocoartf1671
{\fonttbl\f0\fswiss\fcharset0 Helvetica;\f1\fnil\fcharset0 Menlo-Bold;\f2\fnil\fcharset0 Menlo-Regular;
}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;\red251\green2\blue7;\red0\green0\blue0;
\red251\green2\blue7;}
{\*\expandedcolortbl;;\csgray\c0;\cssrgb\c100000\c14913\c0;\cssrgb\c0\c0\c0;
\cssrgb\c100000\c14913\c0;}
\margl1440\margr1440\vieww10800\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 The main program is \'91haplotypePhaser.py\'92. Please execute it like any other python 3 program. For example, in Mac Terminal, after changing to the correct directory, type:\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0

\f1\b \cf2 \CocoaLigature0 python3 haplotypePhaser.py
\f0\b0 \cf0 \CocoaLigature1 \
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0
\cf0 After the user execute the program, the program will print the following prompt: \

\f1\b Please type your file path to start phasing.
\f2\b0 \

\f1\b Input file path (press Enter when finish): 
\f0\b0 \cf3 [please input your path for the genotype data file to read in. E.g.: \'91test_data_1.txt\'92]\

\f1\b \cf0 Output file name (press Enter when finish):
\f2\b0  
\f0 \cf3 [please input your name for the haplotype data file to write out. E.g. \'91test_data_1_sol.txt\'92]\
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0
\cf4 Then the program will read the designated input file, and create and write a file in the same directory, with the designated output file name. The program will keep printing the current running time and percentage of completion (progress) as following:
\f2\fs22 \cf2 \CocoaLigature0 \
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0

\f1\b\fs24 Running Time: 1361.6 s\
Speed:        26.8 SNP/s\
Progress:     67 %
\f0\b0 \cf4 \CocoaLigature1 \
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0
\cf4 The program will end with the following output:\
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f1\b \cf4 OK end!!!\
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\b0 \cf4 \
PS: in my computer, I figured that Terminal runs much slower than other compilers, such as Sublime Text. Therefore, I recommend \cf5 not using Terminal\cf4  because it might cause the program to run over time.}
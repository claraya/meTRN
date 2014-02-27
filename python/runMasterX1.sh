#!/usr/bin/sh

# Generate data quality reports:
python mapData.py --path ~/meTRN/ --mode review --organism ce --source extras --infile cetrn/configure_final_summary.txt --peaks ce_reporting_com_xx_raw --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6,RPC-15:RPC-1,DCC4640_LIN-35:DCC4640_EFL-1,HAM-1:MEP-1,CEH-28:CEH-38
python mapData.py --path ~/meTRN/ --mode revise --organism ce --source extras --infile cetrn/configure_final_summary.txt --peaks ce_reporting_com_xx_raw --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6,RPC-15:RPC-1,DCC4640_LIN-35:DCC4640_EFL-1,HAM-1:MEP-1,CEH-28:CEH-38 --library ce/blacklisted/v2
python mapData.py --path ~/meTRN/ --mode survey --organism ce --source extras --infile cetrn/configure_final_summary.txt --peaks ce_reporting_com_xx_raw --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6,RPC-15:RPC-1,DCC4640_LIN-35:DCC4640_EFL-1,HAM-1:MEP-1,CEH-28:CEH-38
python mapData.py --path ~/meTRN/ --mode report --organism ce --source extras --infile cetrn/configure_final_summary.txt --peaks ce_reporting_com_xx_raw --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6,RPC-15:RPC-1,DCC4640_LIN-35:DCC4640_EFL-1,HAM-1:MEP-1,CEH-28:CEH-38

python mapData.py --path ~/meTRN/ --mode review --organism ce --source extras --infile cetrn/configure_final_summary.txt --peaks ce_finalized_com_cx_raw --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6,RPC-15:RPC-1,DCC4640_LIN-35:DCC4640_EFL-1,HAM-1:MEP-1,CEH-28:CEH-38
python mapData.py --path ~/meTRN/ --mode revise --organism ce --source extras --infile cetrn/configure_final_summary.txt --peaks ce_finalized_com_cx_raw --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6,RPC-15:RPC-1,DCC4640_LIN-35:DCC4640_EFL-1,HAM-1:MEP-1,CEH-28:CEH-38 --library ce/blacklisted/v2
python mapData.py --path ~/meTRN/ --mode survey --organism ce --source extras --infile cetrn/configure_final_summary.txt --peaks ce_finalized_com_cx_raw --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6,RPC-15:RPC-1,DCC4640_LIN-35:DCC4640_EFL-1,HAM-1:MEP-1,CEH-28:CEH-38
python mapData.py --path ~/meTRN/ --mode report --organism ce --source extras --infile cetrn/configure_final_summary.txt --peaks ce_finalized_com_cx_raw --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6,RPC-15:RPC-1,DCC4640_LIN-35:DCC4640_EFL-1,HAM-1:MEP-1,CEH-28:CEH-38

python mapData.py --path ~/meTRN/ --mode review --organism ce --source extras --infile cetrn/configure_final_summary.txt --peaks ce_corrected_com_cx_raw --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6,RPC-15:RPC-1,DCC4640_LIN-35:DCC4640_EFL-1
python mapData.py --path ~/meTRN/ --mode revise --organism ce --source extras --infile cetrn/configure_final_summary.txt --peaks ce_corrected_com_cx_raw --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6,RPC-15:RPC-1,DCC4640_LIN-35:DCC4640_EFL-1
python mapData.py --path ~/meTRN/ --mode survey --organism ce --source extras --infile cetrn/configure_final_summary.txt --peaks ce_corrected_com_cx_raw --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6,RPC-15:RPC-1,DCC4640_LIN-35:DCC4640_EFL-1
python mapData.py --path ~/meTRN/ --mode report --organism ce --source extras --infile cetrn/configure_final_summary.txt --peaks ce_corrected_com_cx_raw --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6,RPC-15:RPC-1,DCC4640_LIN-35:DCC4640_EFL-1

python mapData.py --path ~/meTRN/ --mode review --organism ce --source extras --infile cetrn/configure_final_summary.txt --peaks ce_selection_com_cx_raw --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6,RPC-15:RPC-1
python mapData.py --path ~/meTRN/ --mode revise --organism ce --source extras --infile cetrn/configure_final_summary.txt --peaks ce_selection_com_cx_raw --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6,RPC-15:RPC-1 --library ce/blacklisted/v2
python mapData.py --path ~/meTRN/ --mode survey --organism ce --source extras --infile cetrn/configure_final_summary.txt --peaks ce_selection_com_cx_raw --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6,RPC-15:RPC-1
python mapData.py --path ~/meTRN/ --mode report --organism ce --source extras --infile cetrn/configure_final_summary.txt --peaks ce_selection_com_cx_raw --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6,RPC-15:RPC-1


#top
#bash runMasterX1.sh
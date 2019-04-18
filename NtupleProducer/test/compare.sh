awk '{print $1}' muon_dump.txt > muons_c1ND.txt
awk '{print $1}' muon_IPHC.txt > muons_c1IPHC.txt
diff muons_c1ND.txt muons_c1IPHC.txt > muons_diff.txt

awk '{print $1}' electron_dump.txt > electrons_c1ND.txt
awk '{print $1}' electron_IPHC.txt > electrons_c1IPHC.txt
diff electrons_c1ND.txt electrons_c1IPHC.txt > electrons_diff.txt

awk '{print $1}' tau_dump.txt > taus_c1ND.txt
awk '{print $1}' tau_IPHC.txt > taus_c1IPHC.txt
diff taus_c1ND.txt taus_c1IPHC.txt > taus_diff.txt

awk '{print $1}' jet_dump.txt > jets_c1ND.txt
awk '{print $1}' jet_IPHC.txt > jets_c1IPHC.txt
diff jets_c1ND.txt jets_c1IPHC.txt > jets_diff.txt


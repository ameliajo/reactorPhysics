



import subprocess
#subprocess.check_output(['ls','-l']) #all that is technically needed...
#print(subprocess.check_output(['ls','-l']))
#subprocess.check_output(['ls','-l'])
subprocess.run(['ls'])
#subprocess.run(['mkdir','hello'])
subprocess.run(['python','../GettingXS/realisticCopy/grid_3x3.py'])
subprocess.run(['cp','../GettingXs/realisticCopy/XS.py','../MOC'])
subprocess.run(['python','../MOC/moc.py'])

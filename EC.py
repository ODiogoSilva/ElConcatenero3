#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
#  EC.py
#  
#  Copyright 2013 Unknown <diogo@arch>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

import sys
from PyQt4 import QtGui, QtCore

class EC (QtGui.QMainWindow):

	def __init__(self):
		super(EC, self).__init__()
        
		self.initUI()
        
	def initUI(self): 
    
		## Geometry
		self.setGeometry(300, 300, 290, 150)
		self.center()
		## Header
		self.setWindowTitle('ElConcatenero')   
		#self.setWindowIcon(QtGui.QIcon('sombrero.png')) # ICON MISSING HERE
        
		# Menu bar
		self.createMenus()
        
        
		# Status bar

	def createMenus (self):
		""" Creation of the menus for the menu bar """
		# FILE MENU
		openFile = QtGui.QAction(QtGui.QIcon('open.png'), 'Open', self)
		openFile.setShortcut('Ctrl+O')
		openFile.setStatusTip('Open new File')
		#openFile.triggered.connect(self.showDialog)
        
		menubar = self.menuBar()
		fileMenu = menubar.addMenu('&File')
		fileMenu.addAction(openFile)
		
		# MODE MENU
		
		concatenation = QtGui.QAction(QtGui.QIcon('open.png'),'Concatenation',self)
		concatenation.setShortcut('Ctrl+Alt+C')
		concatenation.setStatusTip('Opens the concatenation window')
		
		modeMenu = menubar.addMenu('&Mode')
		modeMenu.addAction(concatenation)
			
	def center(self):
		""" Function used to center the main window in the Desktop environment """
		qr = self.frameGeometry()
		cp = QtGui.QDesktopWidget().availableGeometry().center()
		qr.moveCenter(cp)
		self.move(qr.topLeft())

	def closeEvent(self, event):
		""" Generates an event at closing prompting the user about closing the program """
		reply = QtGui.QMessageBox.question(self, 'Message',
			"Are you sure to quit?", QtGui.QMessageBox.Yes | 
			QtGui.QMessageBox.No, QtGui.QMessageBox.No)

		if reply == QtGui.QMessageBox.Yes:
			event.accept()
		else:
			event.ignore()

	def keyPressEvent(self, e):
		""" Exit the program using the 'ESC' key """
		if e.key() == QtCore.Qt.Key_Escape:
			self.close()

def main():
    
	app = QtGui.QApplication(sys.argv)
	ex = EC()
	ex.show()
	sys.exit(app.exec_())

if __name__ == '__main__':
	main()

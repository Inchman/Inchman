# Import all TKinter stuff
from Tkinter import Tk, RIGHT, BOTH, RAISED, Text, END
from ttk import Frame, Button, Style, Label
import tkMessageBox, tkFileDialog
import os.path

# The GUI class
class InchmanWindow(Frame):
  
    # constructor
    def __init__(self, parent):
        # call base setting parent window and background
        Frame.__init__(self, parent)   

        # save reference to root window
        self.parent = parent
        
        # create ui
        self.initUI()
    
    # init the user interface
    def initUI(self):        
        # title
        self.parent.title("Inchman process")
        
        # use defualt style
        self.style = Style()
        self.style.theme_use("default")

        # create frame
        frame = Frame(self, relief=RAISED, borderwidth=1)
        frame.pack(fill=BOTH, expand=1)

        # create text box in frame
        self.text = Text(frame)
        self.text.pack(fill=BOTH, expand=1)
        self.pack(fill=BOTH, expand=1)

        # create some tags for differnent colors
        self.text.tag_config("success", foreground="green")

    # adds a text line
    def insert_text(self, text, tag=None):
        self.text.insert(END, text, tag)

class InchmanGUI():
    def __init__(self):
        # create root window
        self.root = Tk()
        self.root.geometry("1000x500+300+300")

        # create the GUI on top of the root window
        self.app = InchmanWindow(self.root)
    
    def getWorkingDirectory(self):
        return tkFileDialog.askdirectory(title="Choose working directory..", parent=self.app, initialdir=os.path.expanduser("~"), mustexist=True) + "/"

    def update(self):
        self.root.update_idletasks()

    def insert_text(self, text, tag=None):
        self.app.insert_text(text, tag)

    def info(self, text):
        tkMessageBox.showinfo("Inchman info", text, parent=self.app)

    def error(self, text):
        tkMessageBox.showerror("Inchman error", text)

import os
import PyQt6 as qt
from PyQt6 import QtGui as qtG
from PyQt6 import QtCore as qtC
from PyQt6 import QtWidgets as qtW
from graphics_items import *
from constantes import * 

class ClickScene (qtW.QGraphicsScene) :
    def __init__ (self, parent_window = None, *args, **kwargs) -> None :
        '''
        init the clickable view
        
        Args :
            parent_window : the clickable view
        '''
        super ().__init__ (*args, **kwargs)
        self.parent_window = parent_window
        
    def mousePressEvent (self, event) -> None :
        '''
        Event triggered when the mouse is pressed in the scene
        
        Args:
            event : the event triggered
        '''
        item = self.itemAt (event.scenePos (), qtG.QTransform ())
        if not isinstance (item, HighlightableRectJob) :
            if self.parent_window :
                self.parent_window.reset_visibility ()
        
        super().mousePressEvent (event)

class ZoomView (qtW.QGraphicsView) :
    def __init__ (self, scene : qtW.QGraphicsScene) -> None :
        '''
        init the zoom view
        
        Args :
            scene (qtW.QGraphicsScene) : the scene to display
        '''
        super ().__init__ (scene)
        self._zoom = 1.0
        self.mode = 0
        self.nb_mode = 2

    def wheelEvent (self, event : qtG.QWheelEvent) -> None :
        '''
        Event triggered when the mouse wheel is scrolled
        
        Args :
            event (qtG.QWheelEvent) : the event triggered
        '''
        if event.modifiers () == qtC.Qt.KeyboardModifier.ControlModifier :
            angle = event.angleDelta ().y ()
            factor = 1.05 if angle > 0 else 1 / 1.05 # speed of the zoom
            self._zoom *= factor
            if self.mode == 0 :
                self.scale (factor, factor)
            elif self.mode == 1 :
                self.scale (factor, 1)
        else :
            super ().wheelEvent (event)

    def reset_zoom (self) -> None :
        '''
        Reset the zoom of the view
        '''
        self.resetTransform ()
        self._zoom = 1.0
    
    def mode_zoom (self, button : qtW.QPushButton) -> None :
        '''
        change the mode of the zoom
        
        Args :
            - button (qtW.QPushButton) : the button that triggered the event
        '''
        if 0 <= self.mode < self.nb_mode - 1 :
            self.mode += 1
        elif self.mode == self.nb_mode - 1 :
            self.mode = 0
          
        if self.mode == 0 :
            image = os.path.join (os.path.dirname (__file__), FILE_ICON, FOUR_DIRECTIONNEL_ARROW_ICON)
            icon = qtG.QIcon (image)
            button.setIcon (icon)
        elif self.mode == 1 :
            image = os.path.join (os.path.dirname (__file__), FILE_ICON, TWO_DIRECTIONNEL_ARROW_ICON)
            icon = qtG.QIcon (image)
            button.setIcon (icon)
import PyQt6 as qt
from PyQt6 import QtGui as qtG
from PyQt6 import QtCore as qtC
from PyQt6 import QtWidgets as qtW
from constantes import *

class CurvedArrow (qtW.QGraphicsPathItem) :
    def __init__ (self, start : qtC.QPointF, end : qtC.QPointF) -> None :
        '''
        init the curve arrow 
        
        Args :
            - start : the start point
            - end : the end point
        '''
        super().__init__ (None)
        self.setAcceptHoverEvents (True)
        self.arrow_size = 10
        self.arrow_height = 6
        self.real_end = end
        self.start = start
        self.end = qtC.QPointF (self.real_end.x () - self.arrow_size, self.real_end.y ())
        self.setZValue (ARROW_ZVALUE)
        self.highlight_pen = qtG.QPen (qtG.QColor (DEFAULT_HIGHTLIGHT_COLOR), 3)
        self.default_pen = qtG.QPen (qtG.QColor (DEFAULT_PEN_COLOR), 2)
        self.brush = qtG.QBrush (qtG.QColor (DEFAULT_BRUSH_COLOR))
        self.setPen (self.default_pen)
        self.setBrush (self.brush)
        self.setPath (self.create_curve_path ())
        
    def create_curve_path (self) -> qtG.QPainterPath :
        '''
        méthode that create the path of the arrow
        
        Return :
            - path (qtG.QPainterPath) : the path of the arrow
        '''
        path = qtG.QPainterPath ()
        path.moveTo (self.start)
        dx = self.real_end.x () - self.start.x ()
        dy = self.real_end.y () - self.start.y ()
        if dx == 0 and dy == 0 :
            val_x = 0
            val_y = 0
        else :
            val_x = 60
            val_y = 0
        ctrl1 = qtC.QPointF (self.start.x () + val_x, self.start.y () + val_y)
        ctrl2 = qtC.QPointF (self.start.x () - val_x, self.end.y () + val_y)
        path.cubicTo (ctrl1, ctrl2, self.end)
        return path
    
    def paint (self, painter : qtG.QPainter, option, widget = None) -> None :
        '''
        méthode that paint the arrow
        
        Args :
            - painter (qtG.QPainter)
            - option 
            - widget
        '''
        painter.setRenderHint (qtG.QPainter.RenderHint.Antialiasing)
        painter.setPen (self.pen ())
        painter.drawPath (self.path ())
        
        tip = self.real_end
        base1 = qtC.QPointF (tip.x () - self.arrow_size, tip.y () - self.arrow_height)
        base2 = qtC.QPointF (tip.x () - self.arrow_size, tip.y () + self.arrow_height)
        triangle = qtG.QPolygonF ([tip, base1, base2])
        painter.setBrush (self.brush)
        painter.drawPolygon (triangle)
        
    def boundingRect (self) -> qtC.QRectF :
        '''
        return the bounding rectangle of the arrow
        
        Return :
            - (qtC.QRectF) : the bounding rectangle of the arrow
        '''
        rect = self.path ().boundingRect ()
        extra = 10
        return rect.adjusted (-extra, -extra, extra, extra)
    
    def hoverEnterEvent (self, event : qtC.QEvent) -> None:
        '''
        Highlight the arrow when the mouse is over
        '''
        self.setPen (self.highlight_pen)
        self.update ()
        super ().hoverEnterEvent (event)
    
    def hoverLeaveEvent (self, event : qtC.QEvent) -> None:
        '''
        Stop highlighting the arrow when the mouse leaves
        '''
        self.setPen (self.default_pen)
        self.update ()
        super ().hoverLeaveEvent (event)

class HighlightableRect (qtW.QGraphicsRectItem) :
    def __init__ (self, *args, parent_window = None, data : dict = None, callback = None, hover_label = None, 
                  default_pen : qtG.QPen = qtG.QPen (qtG.QColor (DEFAULT_PEN_COLOR), 1), highlight_pen : qtG.QPen = qtG.QPen (qtG.QColor (DEFAULT_HIGHTLIGHT_COLOR), 2), clic_pen : qtG.QPen = qtG.QPen (qtG.QColor (DEFAULT_PEN_COLOR), 2, qtC.Qt.PenStyle.DotLine), **kwargs) -> None :
        '''
        Init the highlightable rectangle
        
        Args :
            args : arguments of the rectangle
            parent_window : the clickable view
            data (dict) : the data to display
            callback : function to call when the rectangle is clicked
            hover_label : label to display on hover
            default_pen (qtG.QPen) : the pen of the rectangle
            highlight_pen (qtG.QPen) : the pen of the rectangle when hovered
            kwargs : arguments of the rectangle
        '''
        super ().__init__ (*args, **kwargs)
        self.parent_window = parent_window
        self.data = data
        self.callback = callback
        self.hover_label = hover_label
        self.default_pen = default_pen
        self.save_default_pen = self.default_pen
        self.highlight_pen = highlight_pen
        self.clic_pen = clic_pen
        self.setAcceptHoverEvents (True)
        self.setPen (self.default_pen)
 
    def hoverMoveEvent (self, event : qtC.QEvent) -> None :
        '''
        Event triggered when the mouse is moved on the rectangle
        
        Args :
            event (qtC.QEvent) : the event triggered
        '''
        # Update the position of the label to follow the mouse
        if self.hover_label :
            self.hover_label.move (qtG.QCursor.pos () + qtC.QPoint (10, 10))
            
        super ().hoverMoveEvent (event)

    def hoverLeaveEvent (self, event : qtC.QEvent) -> None :
        '''
        Event triggered when the mouse leaves the rectangle
        
        Args :
            event (qtC.QEvent) : the event triggered
        '''
        self.setPen (self.default_pen)
        if self.hover_label :
            self.hover_label.setVisible (False)
        
        super ().hoverLeaveEvent (event)
        
class HighlightableRectJob (HighlightableRect) :
    def hoverEnterEvent (self, event : qtC.QEvent) -> None :
        '''
        Event triggered when the mouse is on the rectangle
        
        Args :
            event (qtC.QEvent) : the event triggered
        '''
        
        self.setPen (self.highlight_pen)
        
        # Display the job data in a tooltip
        if self.hover_label and self.data :
            info = f"Job : {self.data[NOMJ]}"
            info += f"\nId : {self.data[IDJ]}"
                
            if self.data[IDPARENT] != -1 :
                info += f"\nParent : {self.data[IDPARENT]}"
            
            val = f"{self.data[DEBUT]:,}"
            val = val.replace(",", " ")
            info += f"\nStart : {val}"
            if self.data[FIN] != -1 :
                val = f"{self.data[FIN]:,}"
                val = val.replace(",", " ")
                val1 = f"{self.data[FIN] - self.data[DEBUT]:,}"
                val1 = val1.replace(",", " ")
                info += f"\nEnd : {val}"
                info += f"\nDuration : {val1}"
            if not self.data[SUCCESS] :
                info += f"\nError : {self.data[ERROR]}"
            self.hover_label.setText (info)
            self.hover_label.adjustSize ()
            self.hover_label.move (qtG.QCursor.pos ()) # move the label to the mouse position
            self.hover_label.setVisible (True) # show the label
        
        super ().hoverEnterEvent (event)
    
    def mousePressEvent (self, event : qtC.QEvent) -> None :
        '''
        Event triggered when the mouse is press
        
        Args :
            event (qtC.QEvent) : the event triggered
        '''
        if self.callback and self.data :
            self.callback (self.data[IDJ])
        
        super ().mousePressEvent (event)

class HighlightableRectQuitting (HighlightableRect) :
    def hoverEnterEvent (self, event : qtC.QEvent) -> None :
        '''
        Event triggered when the mouse is on the rectangle
        
        Args :
            event (qtC.QEvent) : the event triggered
        '''
        
        self.setPen (self.highlight_pen)
        
        # Display the quitting data in a tooltip
        if self.hover_label and self.data :
            info = f"Quitting"
            info += f"\nError : {self.data[ERROR]}"
            self.hover_label.setText (info)
            self.hover_label.adjustSize ()
            self.hover_label.move (qtG.QCursor.pos ()) # move the label to the mouse position
            self.hover_label.setVisible (True) # show the label
        
        super ().hoverEnterEvent (event)

class HighlightableRectMutex (HighlightableRect) :
    def hoverEnterEvent (self, event : qtC.QEvent) -> None :
        '''
        Event triggered when the mouse is on the rectangle
        
        Args :
            event (qtC.QEvent) : the event triggered
        '''
        
        self.setPen (self.highlight_pen)
        
        # Display the job data in a tooltip
        if self.hover_label and self.data :
            info = f"Mutex : {self.data[NAMEM]}"
            val = f"{self.data[NEED]:,}"
            val = val.replace(",", " ")
            info += f"\nNeed : {val}"
            if self.data[HAVE] != -1 :
                val = f"{self.data[HAVE]:,}"
                val = val.replace(",", " ")
                info += f"\nHave : {val}"
            if self.data[FREE] != -1 :
                val = f"{self.data[FREE]:,}"
                val = val.replace(",", " ")
                info += f"\nFree : {val}"
        
            self.hover_label.setText (info)
            self.hover_label.adjustSize ()
            self.hover_label.move (qtG.QCursor.pos ()) # move the label to the mouse position
            self.hover_label.setVisible (True) # show the label
        
        super ().hoverEnterEvent (event)

class HighlightableRectTag (HighlightableRect) :
    def hoverEnterEvent (self, event : qtC.QEvent) -> None :
        '''
        Event triggered when the mouse is on the rectangle
        
        Args :
            event (qtC.QEvent) : the event triggered
        '''
        
        self.setPen (self.highlight_pen)
        
        # Display the job data in a tooltip
        if self.hover_label and self.data :
            info = f"{self.data[NAMET]}"
            val = f"{self.data[DT]:,}"
            val = val.replace(",", " ")
            info += f"\nStart : {val}"
            if self.data[FT] != -1 :
                val = f"{self.data[FT]:,}"
                val = val.replace(",", " ")
                val1 = f"{self.data[FT] - self.data[DT]:,}"
                val1 = val1.replace(",", " ")
                info += f"\nStop : {val}"
                info += f"\nDuration : {val1}"
        
            self.hover_label.setText (info)
            self.hover_label.adjustSize ()
            self.hover_label.move (qtG.QCursor.pos ()) # move the label to the mouse position
            self.hover_label.setVisible (True) # show the label
        
        super ().hoverEnterEvent (event)
        
class HighlightableRectDeadlock (HighlightableRect) :
    def hoverEnterEvent (self, event : qtC.QEvent) -> None :
        '''
        Event triggered when the mouse is on the rectangle
        
        Args :
            event (qtC.QEvent) : the event triggered
        '''
        
        self.setPen (self.highlight_pen)
        
        # Display the deadlock data
        if self.hover_label and self.data :
            deadlock = self.data[DEADLOCK]
            threads = self.data[THREADS]
            val = f"{deadlock[DEADLOCK_TIME]:,}"
            val = val.replace(",", " ")
            info = f"⚠️ DEADLOCK SITUATION at {val} ms ⚠️"
            for th in deadlock[DEADLOCK_CYCLE] :
                for thread in threads :
                    if thread[IDT] == th[DEADLOCK_TID] :
                        thread_name = thread[NOMT]
                info += f"\n        The thread \"{thread_name}\""
                have = ""
                for mutex in th[DEADLOCK_HAVE] :
                    have += mutex + " "
                info += f"\n            have : {have}"
                info += f"\n            need : {th[DEADLOCK_NEED]}"
        
            self.hover_label.setText (info)
            self.hover_label.adjustSize ()
            self.hover_label.move (qtG.QCursor.pos ()) # move the label to the mouse position
            self.hover_label.setVisible (True) # show the label
        
        super ().hoverEnterEvent (event)

class TriangleSlider (qtW.QGraphicsPolygonItem) :
    def __init__ (self, left_margin : int, top_margin : int, pixels_per_unit : int, thread_h : int, nb_threads : int, slider : qtW.QSlider, x_min : int, x_max : int) -> None :
        '''
        Init the triangle slider
        
        Args :
            left_margin (int) : the left margin of the view
            top_margin (int) : the top margin of the view
            pixels_per_unit (int) : the number of pixels per unit
            thread_h (int) : the height of the threads
            nb_threads (int) : the number of threads
            slider (qtW.QSlider) : the slider to move
            x_min (int) : the minimum value of the slider
            x_max (int) : the maximum value of the slider
        '''
        points = [qtC.QPointF (left_margin + pixels_per_unit, top_margin), 
                  qtC.QPointF (left_margin + pixels_per_unit - 5, top_margin - 10), 
                  qtC.QPointF (left_margin + pixels_per_unit + 5, top_margin - 10)]
        super ().__init__ (qtG.QPolygonF (points))
        self.left_margin = left_margin
        self.top_margin = top_margin
        self.pixels_per_unit = pixels_per_unit
        self.thread_h = thread_h
        self.nb_threads = nb_threads
        self.slider = slider
        self.x_min = x_min
        self.x_max = x_max * self.pixels_per_unit
        
        self.setBrush (qtG.QColor (DEFAULT_BRUSH_COLOR))
        self.setZValue (TRIANGLE_SLIDER_ZVALUE)
        self.setFlag (qtW.QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        self.setFlag (qtW.QGraphicsItem.GraphicsItemFlag.ItemSendsScenePositionChanges)
        self.setCursor (qtC.Qt.CursorShape.OpenHandCursor)

    def mouseMoveEvent (self, event : qtC.QEvent) -> None :
        '''
        Event triggered when the mouse move the slider triangle
        
        Args : 
            event (qtC.QEvent) : the event triggered
        '''
        new_x = event.scenePos ().x () - self.left_margin # get the x position of the mouse
        
        clamped_x = int (max (self.x_min, min (new_x, self.x_max)) / self.pixels_per_unit) * self.pixels_per_unit
        
        self.slider.setValue (clamped_x // self.pixels_per_unit) # update the slider value => this will trigger the update_slider method

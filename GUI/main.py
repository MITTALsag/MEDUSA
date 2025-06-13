import os
import sys
import json
import PyQt6 as qt
from PyQt6 import QtGui as qtG
from PyQt6 import QtCore as qtC
from PyQt6 import QtWidgets as qtW
from graphics_items import *
from graphics_view import *
from constantes import * 

################################################################
# -- Functions to create generic buttons in the GUI --
def create_button (text : str, callback, checkable : bool = False) -> qtW.QPushButton :
    '''
    Helper to create text buttons
    
    Args :
        text (str) : the text of the button
        callback : the function associate to the button
        checkable (bool) : if the button is checkable or not (default is False)
        
    Return :
        button (qtW.QPushButton) : the button create
    '''
    button = qtW.QPushButton (text)
    if checkable :
        button.setCheckable (True)
        button.setChecked (False)
    button.clicked.connect (callback)
    
    background_color = "#{:06X}".format (BUTTON_BACKGROUND_COLOR)
    border_color = "#{:06X}".format (BUTTON_BORDER_COLOR)
    hover_color = "#{:06X}".format (BUTTON_HOVER_COLOR)
    checked_color = "#{:06X}".format (BUTTON_CHECKED_COLOR)
    hover_text_color = "#{:06X}".format (BUTTON_HOVER_TEXT_COLOR)
    checked_text_color = "#{:06X}".format (BUTTON_CHECKED_TEXT_COLOR)
    button.setStyleSheet (f'''
        QPushButton {{
            background-color : {background_color};
            border : 1px solid {border_color};
            padding : 5px;
            border-radius : 10px;
        }}
        QPushButton:checked {{
            background-color : {checked_color};
            color : {checked_text_color};
        }}
        QPushButton:hover {{
            background-color : {hover_color};
            color : {hover_text_color};
        }}
    ''')
    button.setCursor (qtC.Qt.CursorShape.PointingHandCursor)
    return button

def create_image_button (image_path : str, callback, checkable : bool = False) -> qtW.QPushButton :
    '''
    Helper to create icon buttons
    
    Args :
        image_path (str): the path to the image file
        callback : the function associate to the button
        checkable (bool) : if the button is checkable or not (default is False)
        
    Return :
        button (qtW.QPushButton) : the button create
    '''
    button = qtW.QPushButton ()
    if checkable :
        button.setCheckable (True)
        button.setChecked (False)
    button.clicked.connect (callback)

    # Set the icon for the button
    image = os.path.join (os.path.dirname (__file__), FILE_ICON, image_path)
    icon = qtG.QIcon (image)
    button.setIcon (icon)

    # Optionally, set the size of the icon
    button.setIconSize (button.sizeHint ())
    button.setFixedSize (button.iconSize () + qtC.QSize (0, 6))
    
    background_color = "#{:06X}".format (BUTTON_BACKGROUND_COLOR)
    hover_color = "#{:06X}".format (BUTTON_HOVER_COLOR)
    checked_color = "#{:06X}".format (BUTTON_CHECKED_COLOR)
    button.setStyleSheet (f'''
        QPushButton {{
            background-color : {background_color};
            border-radius : 2px;
        }}
        QPushButton:checked {{
            background-color : {checked_color};
        }}
        QPushButton:hover {{
            background-color : {hover_color};
        }}
    ''')
    button.setCursor (qtC.Qt.CursorShape.PointingHandCursor)
    return button

def increment_button (title : str, init_value : int, value_min : int, value_max : int, minus_callback, plus_callback, update_callback) -> (qtW.QWidget, qtW.QLineEdit) :
    '''
    Create the button that can change a value
    
    Args :
        - title (str) : the title of the button
        - init_value (int) : the initial value of the button
        - value_min (int) : the minimum value of the button
        - value_max (int) : the maximum value of the button
        - minus_callback : the function to call when the minus button is clicked
        - plus_callback : the function to call when the plus button is clicked
        - update_callback : the function to call when the value is changed in the text area
    
    Returns :
        - button (qtW.QWidget) : the widget that contains the buttons and the value
        - text_value (qtW.QLineEdit) : the text area that contains the value
    '''
    button = qtW.QWidget ()
    button.setLayout (qtW.QGridLayout ())
    button.layout ().setSpacing (0)
    
    background_color = "#{:06X}".format (BUTTON_BACKGROUND_COLOR)
    border_color = "#{:06X}".format (BUTTON_BORDER_COLOR)
    button.setStyleSheet (f'''
        QWidget {{
            background-color : {background_color};
            border : 1px solid {border_color};
            padding : 0px;
            border-radius : 10px;
            margin : 0px;
        }}
        QPushButton {{
            border : none;
            padding : 0px 0px;
        }}
    ''')
    
    # Label for title
    label = qtW.QLabel (title)
    label.setAlignment (qtC.Qt.AlignmentFlag.AlignCenter)
    label.setStyleSheet ("border : 0px;")
    button.layout ().addWidget (label, 0, 0, 1, 3)
    
    # Decrease button
    decrease_button = create_button ("-", minus_callback)
    
    hover_text_color = "#{:06X}".format (BUTTON_HOVER_TEXT_COLOR)
    decrease_button.setStyleSheet (f'''
        QPushButton {{
            background-color : {TRANSPARENT_COLOR};
            border : 0px;
            padding : 0px 0.5px 3.5px 1.5px;
            font-size : 25px;
            border-radius : 10px;
            width : 20px;
            height : 20px;
        }}
        QPushButton:hover {{
            color : {hover_text_color};
        }}
    ''')
    button.layout ().addWidget (decrease_button, 1, 0, 1, 1)
    
    # Text area
    text_value = qtW.QLineEdit (str (init_value))
    border_color = "#{:06X}".format (BUTTON_BORDER_COLOR)
    text_value.setStyleSheet (f"border : 1px solid {border_color};")
    text_value.setFixedWidth (50)
    text_value.setValidator (qtG.QIntValidator (value_min, value_max))
    text_value.editingFinished.connect (update_callback)
    text_value.setAlignment (qtC.Qt.AlignmentFlag.AlignCenter)  # Center the text inside the QLineEdit
    button.layout ().addWidget (text_value, 1, 1, 1, 1)

    # Increase button
    increase_button = create_button ("+", plus_callback)
    
    hover_text_color = "#{:06X}".format (BUTTON_HOVER_TEXT_COLOR)
    increase_button.setStyleSheet (f'''
        QPushButton {{
            background-color : {TRANSPARENT_COLOR};
            border : 0px;
            padding : 0px 0.5px 3.5px 1.5px;
            font-size : 25px;
            border-radius : 10px;
            width : 20px;
            height : 20px;
        }}
        QPushButton:hover {{
            color : {hover_text_color};
        }}
    ''')
    button.layout ().addWidget (increase_button, 1, 2, 1, 1)
    
    return button, text_value
################################################################

################################################################
# -- Create a copy of a HighlightableRect element for zoom effect --
def create_zoomed_copy (original) -> HighlightableRect :
    '''
    Create a copy of the original HighlightableRect element

    Args :
        - original : the HighlightableRect element to copy

    Returns :
        - copy (HighlightableRect) : the zoomed copy of the element
    '''
    # Get the rectangle of the original element
    rect = original.rect ()

    # Dictionary mapping class to its constructor
    rect_classes = {
        HighlightableRectJob : HighlightableRectJob,
        HighlightableRectQuitting : HighlightableRectQuitting,
        HighlightableRectMutex : HighlightableRectMutex,
        HighlightableRectTag : HighlightableRectTag,
        HighlightableRectDeadlock : HighlightableRectDeadlock
    }

    ctor = None
    # Find the constructor corresponding to the original's class
    for cls , constructor in rect_classes.items () :
        if isinstance (original , cls) :
            ctor = constructor
            break

    if ctor is None :
        raise TypeError (f"Unsupported type for zoomed copy : {type (original)}")

    # Prepare arguments for the constructor
    kwargs = {
        "parent_window" : getattr (original, "parent_window", None),
        "data" : getattr (original, "data", None) ,
        "hover_label" : getattr (original, "hover_label", None),
        "default_pen" : original.pen (),
        "highlight_pen" : getattr (original, "highlight_pen", None),
    }

    # Add the callback if present (for Job and Quitting)
    if hasattr (original, "callback") :
        kwargs["callback"] = original.callback

    # Create the copy with the same dimensions and properties
    copy = ctor (rect.x (), rect.y (), rect.width (), rect.height (), **kwargs)
    copy.setBrush (original.brush ())
    copy.setPen (original.pen ())
    copy.setZValue (original.zValue () + 10)  # Ensure it is above the original

    return copy
################################################################

class MainWindow (qtW.QMainWindow) :
    def __init__ (self) -> None :
        '''
        Initialize the main window of the application
        '''
        super ().__init__ ()
        # -- Initialize the atributes of the class --
        self.main_layout = None
        self.menu_layout = None
        self.pos_menu = 0
        self.data_layout = None
        
        # -- Initialize the GUI -- 
        self.setup_main_window ()
        self.setup_menu ()
        self.setup_data_display ()
        
        # listes of widget add on the menu and on the data
        self.list_menu_widgets = []
        self.list_data_widgets = []
    
    ################################################################
    # -- Setup window title and main layout --
    def setup_main_window (self) -> None :
        '''
        Initialize main window settings
        '''   
        self.setWindowTitle ("Chaos project")
        self.showMaximized ()
        main_widget = qtW.QWidget ()
        self.setCentralWidget (main_widget)
        self.main_layout = qtW.QGridLayout (main_widget)
        self.main_layout.setContentsMargins (0, 0, 0, 0)
        self.main_layout.setSpacing (0)
    
    # -- Setup menu layout --
    def setup_menu (self) -> None :
        '''
        Initialize menu section
        '''
        menu = qtW.QWidget ()
        color_hex = "#{:06X}".format (MENU_BACKGROUND_COLOR)
        menu.setStyleSheet (f"background-color: {color_hex};")
        self.menu_layout = qtW.QGridLayout (menu)
        self.menu_layout.setContentsMargins (5, 10, 5, 10)
        self.menu_layout.setSpacing (15)
        self.main_layout.addWidget (menu, 0, 0, 2, 1) # (widget, row, column, row span, column span)
        
        # Menu label
        menu_label = qtW.QLabel ("Menu")
        self.menu_layout.addWidget (menu_label, self.pos_menu, 0, 1, 2, qtC.Qt.AlignmentFlag.AlignCenter)
        self.pos_menu += 1
        
        # Load JSON button
        load_button = create_button ("Load file", self.load_json)
        self.menu_layout.addWidget (load_button, self.pos_menu, 0, 1, 2)
        self.pos_menu += 1

        # Spacer
        spacer_menu = qtW.QSpacerItem (1, 0, qtW.QSizePolicy.Policy.Minimum, qtW.QSizePolicy.Policy.Expanding)
        self.menu_layout.addItem (spacer_menu, 100000, 0)
    
    # -- Setup data layout --
    def setup_data_display (self) -> None :
        '''
        Initialize data display area
        '''
        data_W = qtW.QWidget ()
        color_hex = "#{:06X}".format (DATA_BACKGROUND_COLOR)
        data_W.setStyleSheet (f"background-color: {color_hex};")
        self.data_layout = qtW.QGridLayout (data_W)
        self.data_layout.setContentsMargins (0, 0, 0, 0)
        self.data_layout.setSpacing (0)
        self.main_layout.addWidget (data_W, 0, 1, 2, 10) # (widget, row, column, row span, column span)
    ################################################################
    
    ################################################################
    # -- Getters methods --
    def get_main_layout (self) -> qtW.QGridLayout :
        '''
        Get the main layout of the window
        
        Returns :
            main_layout (qtW.QGridLayout) : the main layout of the window
        '''
        return self.main_layout
    
    def get_menu_layout (self) -> (qtW.QGridLayout, int) :
        '''
        Get the menu layout of the window
        
        Returns :
            menu_layout (qtW.QGridLayout) : the menu layout of the window
            pos_menu (int) : the first free position of the menu in the layout
        '''
        return self.menu_layout, self.pos_menu
    
    def get_data_layout (self) -> qtW.QGridLayout :
        '''
        Get the data layout of the window
        
        Returns :
            data_layout (qtW.QGridLayout) : the data layout of the window
        '''
        return self.data_layout
    ################################################################
    
    ################################################################
    # -- Add methodes --
    def add_button_to_menu_layout (self, button : qtW.QWidget) -> None :
        '''
        Add a button to the menu layout
        
        Args :
            - button (qtW.QWidget) : the button to add to the menu
        '''
        self.menu_layout.addWidget (button, self.pos_menu, 0, 1, 2)
        self.list_menu_widgets.append (button)
        self.pos_menu += 1
        
    def add_widget_to_data_layout (self, widget : qtW.QWidget, row : int, column : int, row_span : int = 1, column_span : int = 1) -> None :
        '''
        Add a widget to the data layout
        
        Args :
            - widget (qtW.QWidget) : the widget to add to the data layout
            - row (int) : the row where to add the widget
            - column (int) : the column where to add the widget
            - row_span (int) : the number of rows that the widget will span (default is 1)
            - column_span (int) : the number of columns that the widget will span (default is 1)
        '''
        self.data_layout.addWidget (widget, row, column, row_span, column_span)
        self.list_data_widgets.append (widget)
    ################################################################
    
    ################################################################
    # -- Load a json file and display it in the GUI --
    def load_json (self) -> None :
        '''
        Load a json file and display it in the GUI
        '''
        file_name, _ = qtW.QFileDialog.getOpenFileName (self, "Charger fichier JSON", "", "JSON Files (*.json)") # Open a file dialog to select a json file
        if file_name : # if a file is selected 
            for widget in self.list_menu_widgets :
                self.menu_layout.removeWidget (widget)
                widget.setParent (None)
                self.pos_menu -= 1
            for widget in self.list_data_widgets :
                self.data_layout.removeWidget (widget)
                widget.setParent (None)
            self.list_menu_widgets.clear ()
            self.list_data_widgets.clear ()
            self.setWindowTitle (file_name)
            with open (file_name, 'r') as f : # open the json file in read mode
                data = json.load (f) # load the json data
                DispalyJsonData (self, data) # display the data in the GUI
    ################################################################

class DispalyJsonData () :
    def __init__ (self, window : MainWindow, data : dict) -> None :
        '''
        Display the json data in the window and add all the objects we need to use the data
        
        Args : 
            - window (MainWindow) : the main window we need to display the data
            - data (dict) : the json data to display
        '''
        self.window = window
        self.data = data
        
        self.final_time = data[FINAL_TIME]
        self.nb_threads = data[NB_THREADS]
        self.nb_jobs_OK = data[NB_JOBS_OK]
        self.nb_jobs_KO = data[NB_JOBS_KO]
        self.nb_jobs_running = data[NB_JOBS_RUNNING]
        self.nb_quitting = data[NB_QUITTING]

        self.threads = data[THREADS] # dictionary with the threads informations
        
        self.deadlock = data[DEADLOCK] # if their is a deadlock situation, the informations are here else this is None

        self.top_margin = TOP_MARGIN_VALUE
        self.left_margin = LEFT_MARGIN_VALUE
        self.pixels_per_unit = PIXELS_PER_UNIT_VALUE
        
        screen_height = qtW.QApplication.primaryScreen ().size ().height ()
        self.thread_h = max (THREADS_HEIGHT_MIN_VALUE, (screen_height - self.top_margin) / (self.nb_threads + 4))
        
        # Create the liste of the item
        self.list_job_item = [None] * (self.nb_jobs_OK + self.nb_jobs_KO + self.nb_jobs_running) # one item := rectangle of the job, item at id i is at pos i
        self.list_quitting_item = [] # one item := rectage of quitting
        self.list_mutex_item = [] # one item := (rectangle of need to have, rectangle of have to free) rectangle have to free can be None if there is no have time
        self.list_tag_item = [] # one item := rectangle of tag
        self.list_deadlock_item = [] # one item := rectagle of deadlock
        self.list_arrow = [] # one := a arrow between two jobs who have a parent link
        
        # Create the tree of the jobs 
        self.list_job_parent = [[None, -1, []] for i in range (self.nb_jobs_OK + self.nb_jobs_KO + self.nb_jobs_running)] # one item := [thread id, parent id, list of children ids]
        
        # Current job
        self.current_job_id = None
        
        # Init the scene and the view
        self.scene = None
        self.view = None
        self.slider = None
        self.slider_t = None
        self.info_widget = None
        self.tracking_rect = None
        self.track_button = None
        self.width_track_rectangle_button = None
        self.width_track_rectangle = INIT_WIDTH_TRACK_RECTAGLE_VALUE
        self.width_track_rectangle_slider = None
        self.width_track_rectangle_slider_label = None
        self.track_mask = None
        self.original_visibility = {}
        self.zoom_factor_button = None
        self.zoom_factor = INIT_ZOOM_VALUE
        self.zoom_slider = None
        self.zoom_slider_label = None
        graphic_area = self.create_graphics_area ()
        
        splitter = qtW.QSplitter (qtC.Qt.Orientation.Vertical) # Create a splitter object
        splitter.addWidget (graphic_area)
        splitter.addWidget (self.info_widget)
        splitter.setSizes ([900, 100])
        self.window.add_widget_to_data_layout (splitter, 0, 0, 1, 10)
        
        # Add the buttons in the menu
        self.mutex_button = None
        self.mutex_colors = {}
        self.arrow_button = None
        self.nb_parent = INIT_PARENT_PATH_VALUE
        self.nb_children = INIT_CHILD_PATH_VALUE
        self.text_value_parent = None
        self.text_value_children = None
        self.add_buttons_to_menu ()
        
        # Init the hover label
        self.hover_label = None
        self.create_hover_label ()
        
        # Display the data
        self.display_time_line ()
        self.display_threads_row ()
        self.display_threads ()
        if self.deadlock is not None :
            self.display_deadlock ()
            
        # Active the loupe
        self.track_button.setChecked (True)
        self.toggle_tracking_rectangle (True)
    
    ################################################################
    # -- Create the view and the scene --
    def create_graphics_area (self) -> qtW.QWidget :
        '''
        Create the graphics area
        
        Returns :
            - graphic_area (qtW.Widget) : the widget that contain all the graphic area
        '''
        # Create the title with the legende
        title = qtW.QLabel (f''' Number of threads : {self.nb_threads} {(f'| Jobs OK : {self.nb_jobs_OK} {JOB_OK_EMOJI} ') if self.nb_jobs_OK > 0 else ''}{(f'| Jobs KO : {self.nb_jobs_KO} {JOB_KO_EMOJI} ') if self.nb_jobs_KO > 0 else ''}{(f'| Jobs unfinished : {self.nb_jobs_running} {JOB_RUNNING_EMOJI} ') if self.nb_jobs_running > 0 else ''}| Mutex {MUTEX_EMOJI} | Tag {TAG_EMOJI} {(f'| Force quit {QUITTING_EMOJI} ') if self.nb_quitting > 0 else ''}{(f'| Deadlock {DEADLOCK_EMOJI} ') if self.deadlock is not None else ''}''')
        
        # Spacer
        spacer_title = qtW.QSpacerItem (1, 0, qtW.QSizePolicy.Policy.Expanding, qtW.QSizePolicy.Policy.Minimum)
        
        # Create the scene and the view
        width = 2 * self.left_margin + (self.final_time + 2) * self.pixels_per_unit
        height = 1.1 * self.top_margin + (self.nb_threads + 0.5) * self.thread_h

        self.scene = ClickScene (parent_window = self)
        self.view = ZoomView (self.scene)
        self.scene.setSceneRect (0, 0, width, height) # set the scene rect to the size of the view (x, y, width, height)
        
        # Change mode zoom button
        mode_zoom_button = create_image_button (FOUR_DIRECTIONNEL_ARROW_ICON, lambda : self.view.mode_zoom (mode_zoom_button))
        
        # Reset zoom button
        reset_zoom_button = create_image_button (RESET_ZOOM_ICON, lambda : self.view.reset_zoom ())
        
        # Create the slider
        self.create_slider ()
        
        # Create the tracking rectangle and the button associeted
        self.create_tracking_rectangle ()
        self.zoom_factor_button = self.create_zoom_factor_button ()
        self.width_track_rectangle_button = self.create_width_track_rectangle_button ()
        self.track_button = create_image_button (GENIUS_EFFECT_ICON, self.toggle_tracking_rectangle, True)
        
        # widget which containe all the graphics items : the title, the loupe button and the zoom button and the slider
        graphic_area = qtW.QWidget ()
        graphic_area.setLayout (qtW.QGridLayout ())
        graphic_area.layout ().setContentsMargins (0, 0, 0, 0)
        graphic_area.layout ().setSpacing (0)
        
        # Add the title and buttons widgets
        graphic_area.layout ().addWidget (title, 0, 0, 1, 1)
        graphic_area.layout ().addItem (spacer_title, 0, 2, 1, 1)
        graphic_area.layout ().addWidget (self.width_track_rectangle_button, 0, 5, 1, 1)
        graphic_area.layout ().addWidget (self.zoom_factor_button, 0, 6, 1, 1)
        graphic_area.layout ().addWidget (self.track_button, 0, 7, 1, 1)
        graphic_area.layout ().addWidget (mode_zoom_button, 0, 8, 1, 1)
        graphic_area.layout ().addWidget (reset_zoom_button, 0, 9, 1, 1)
        
        # Add the other widgets
        graphic_area.layout ().addWidget (self.view, 1, 0, 1, 10)
        graphic_area.layout ().addWidget (self.slider, 2, 0, 1, 10)
        
        return graphic_area
        
    # -- Add the buttons in the menu bar --
    def add_buttons_to_menu (self) -> None :
        '''
        Add the buttons to the menu bar
        '''
        # Mutex button
        self.mutex_button = create_button ("Mutex view", self.toggle_mutex_view, True)
        self.mutex_button.setChecked (False)
        self.window.add_button_to_menu_layout (self.mutex_button)
        
        # Button to enable/disable arrows
        self.arrow_button = create_button ("Heredity view", self.toggle_arrows, True)
        self.arrow_button.setChecked (False)
        self.window.add_button_to_menu_layout (self.arrow_button)
        
        # Parent child button
        parent_child_button = self.create_parent_child_button ()
        self.window.add_button_to_menu_layout (parent_child_button)
    
    # -- Create the parent and child button --
    def create_parent_child_button (self) -> qtW.QWidget :
        '''
        Create the parent child button
        
        Returns :
            - button (qtW.QWidget) : the widget which represents the parent and child button
        '''
        # Parent show button
        parent_button, self.text_value_parent = increment_button ("Number of\nparents", 
                                                                    self.nb_parent, 
                                                                    MIN_PATH_VALUE, 
                                                                    MAX_PATH_VALUE, 
                                                                    lambda : self.decrease_parent_show (MIN_PATH_VALUE), 
                                                                    lambda : self.increase_parent_show (MAX_PATH_VALUE), 
                                                                    self.update_parent_show)
        parent_button.setStyleSheet ('border : 0px')
        
        # Child show button
        child_button, self.text_value_children = increment_button ("Number of\nchildren", 
                                                                    self.nb_children, 
                                                                    MIN_PATH_VALUE, 
                                                                    MAX_PATH_VALUE, 
                                                                    lambda : self.decrease_child_show (MIN_PATH_VALUE), 
                                                                    lambda : self.increase_child_show (MAX_PATH_VALUE), 
                                                                    self.update_child_show)
        child_button.setStyleSheet ('border : 0px')
        
        # Add the child and parent button
        background_color = "#{:06X}".format (BUTTON_BACKGROUND_COLOR)
        border_color = "#{:06X}".format (BUTTON_BORDER_COLOR)
        button = qtW.QWidget ()
        button.setStyleSheet (f'''
            QWidget {{
                background-color : {background_color};
                border : 1px solid {border_color};
                padding : 0px;
                border-radius : 10px;
                margin : 0px;
            }}
        ''')
        button.setLayout (qtW.QGridLayout ())
        button.layout ().setSpacing (0)
        button.layout ().setContentsMargins (1, 1, 1, 1)
        button.layout ().addWidget (child_button, 0, 0)
        button.layout ().addWidget (parent_button, 1, 0)
        
        return button
    
    # -- Create width_track_rectangle_button --
    def create_width_track_rectangle_button (self) -> qtW.QWidget :
        '''
        Create the width track rectangle slider button
        
        Returns :
            - button (qtW.QWidget) : the widget which represents the width track rectangle slider button
        '''
        self.width_track_rectangle_slider = qtW.QSlider (qtC.Qt.Orientation.Horizontal)
        self.width_track_rectangle_slider.setMinimum (MIN_ZOOM_VALUE)
        self.width_track_rectangle_slider.setMaximum (MAX_ZOOM_VALUE)
        self.width_track_rectangle_slider.setMinimumWidth (100)
        self.width_track_rectangle_slider.setMaximumWidth (100)
        self.width_track_rectangle_slider.setMinimumHeight (24)
        self.width_track_rectangle_slider.setMaximumHeight (24)
        self.width_track_rectangle_slider.setValue (self.width_track_rectangle)
        self.width_track_rectangle_slider.valueChanged.connect (self.update_width_track_rectangle)
        
        self.width_track_rectangle_slider_label = qtW.QLabel (f"Range")
        
        widget = qtW.QWidget ()
        layout = qtW.QVBoxLayout (widget)
        layout.setSpacing (0)
        layout.setContentsMargins (0, 0, 0, 0)
        layout.addWidget (self.width_track_rectangle_slider_label, 0, qtC.Qt.AlignmentFlag.AlignCenter)
        layout.addWidget (self.width_track_rectangle_slider, 0, qtC.Qt.AlignmentFlag.AlignCenter)
        
        widget.setVisible (False)
        
        return widget
    
    # -- Create zoom_factor_button --
    def create_zoom_factor_button (self) -> qtW.QWidget :
        '''
        Create the zoom factor slider button
        
        Returns :
            - widget (qtW.QWidget) : the widget which represents the zoom factor slider button
        '''
        self.zoom_slider = qtW.QSlider (qtC.Qt.Orientation.Horizontal)
        self.zoom_slider.setMinimum (MIN_ZOOM_VALUE)
        self.zoom_slider.setMaximum (MAX_ZOOM_VALUE)
        self.zoom_slider.setMinimumWidth (100)
        self.zoom_slider.setMaximumWidth (100)
        self.zoom_slider.setMinimumHeight (24)
        self.zoom_slider.setMaximumHeight (24)
        self.zoom_slider.setValue (self.zoom_factor)
        self.zoom_slider.valueChanged.connect (self.update_zoom_factor)
        
        self.zoom_slider_label = qtW.QLabel (f"Zoom : {self.zoom_factor}x")
        
        widget = qtW.QWidget ()
        layout = qtW.QVBoxLayout (widget)
        layout.setSpacing (0)
        layout.setContentsMargins (0, 0, 0, 0)
        layout.addWidget (self.zoom_slider_label, 0, qtC.Qt.AlignmentFlag.AlignCenter)
        layout.addWidget (self.zoom_slider, 0, qtC.Qt.AlignmentFlag.AlignCenter)
        
        widget.setVisible (False)
        
        return widget
    
    # -- Create the hover label --
    def create_hover_label (self) -> None :
        '''
        Create the hover label to display the data
        '''
        background_color = "#{:06X}".format (TOOLTIP_BACKGROUND_COLOR)
        border_color = "#{:06X}".format (TOOLTIP_BORDER_COLOR)
        
        self.hover_label = qtW.QLabel ("", self.window)
        self.hover_label.setStyleSheet (f"background-color : {background_color}; border : 1px solid {border_color}; padding : 5px; font-size : 10px;")
        self.hover_label.setWindowFlags (qtC.Qt.WindowType.ToolTip) # make the label a tooltip
        self.hover_label.setVisible (False)
    ################################################################
    
    ################################################################
    # -- Methods to create and manege the simple slider --
    def create_slider (self) -> None :
        '''
        create the slider
        '''
        slider_line = qtW.QGraphicsLineItem (self.left_margin + self.pixels_per_unit, 
                                            self.top_margin, 
                                            self.left_margin + self.pixels_per_unit, 
                                            self.top_margin + (self.nb_threads + 0.5) * self.thread_h)
        slider_line.setPen (qtG.QPen (qtG.QColor (BLACK_COLOR), 1, qtC.Qt.PenStyle.DotLine))
        slider_line.setZValue (SLIDER_LINE_ZVALUE)
        self.scene.addItem (slider_line)
        
        time_label = qtW.QGraphicsSimpleTextItem ("0")
        time_label.setPos (self.left_margin + self.pixels_per_unit - time_label.boundingRect ().width () / 2, self.top_margin - 25)
        self.scene.addItem (time_label)
        
        self.slider = qtW.QSlider (qtC.Qt.Orientation.Horizontal)
        self.slider.setMinimum (0)
        self.slider.setMaximum (self.final_time)
        self.slider.valueChanged.connect (lambda value : self.update_slider (value, slider_line, time_label))
        
        self.slider_t = TriangleSlider (self.left_margin, self.top_margin, self.pixels_per_unit, self.thread_h, self.nb_threads, self.slider, 0, self.final_time)
        self.scene.addItem (self.slider_t)
        
        self.info_widget = qtW.QTextEdit ()
        self.info_widget.setReadOnly (True)  # Make the text widget read-only
        
        self.update_slider (0, slider_line, time_label) # update the slider to the first value
    
    def update_slider(self, value : int, slider_line : qtW.QGraphicsLineItem, time_label : qtW.QGraphicsSimpleTextItem) -> None :
        '''
        Update the position of the slider and check for intersections with jobs
        
        Args:
            - value (int): The value of the slider
            - slider_line (qtW.QGraphicsLineItem) : the line that slide on the scene
            - time_label (qtW.QGraphicsSimpleTextItem) : the label that slide on the scene
        '''
        center_x = self.left_margin + (value + 1) * self.pixels_per_unit

        # Update slider position
        self.slider_t.setPos (value * self.pixels_per_unit, 0)
        slider_line.setLine (center_x, self.top_margin, center_x, self.top_margin + (self.nb_threads + 0.5) * self.thread_h)
        slider_line.setZValue (SLIDER_LINE_ZVALUE)
        val = f"{value:,}"
        val = val.replace(",", " ")
        time_label.setText (f"{val}")
        time_label.setPos (self.left_margin + (value + 1) * self.pixels_per_unit - time_label.boundingRect ().width () / 2, self.top_margin - 25)
        self.display_active_info (value)

        # Update the tracking rectangle position
        if self.tracking_rect is not None :
            self.update_tracking_rectangle (value)
    
    # -- Methods to display the active information at the current time --
    def display_active_info (self, value : int) -> None :
        '''
        Display the active info at the current time
        
        Args :
            value (int) : the current time
        '''
        self.info_widget.clear ()
        val = f"{value:,}"
        val = val.replace(",", " ")
        self.info_widget.append (f"🕒 Time : {val} ms\n")
        
        if self.deadlock is not None :
            val = f"{self.deadlock[DEADLOCK_TIME]:,}"
            val = val.replace(",", " ")
            self.info_widget.append (f"⚠️ <u><b>DEADLOCK SITUATION at {val} ms</b></u> ⚠️")
            for th in self.deadlock[DEADLOCK_CYCLE] :
                for thread in self.threads :
                    if thread[IDT] == th[DEADLOCK_TID] :
                        thread_name = thread[NOMT]
                self.info_widget.append (f"        The thread \"{thread_name}\"")
                have = ""
                for mutex in th[DEADLOCK_HAVE] :
                    have += mutex + " "
                self.info_widget.append (f"            have : {have}")
                self.info_widget.append (f"            need : {th[DEADLOCK_NEED]}")
            self.info_widget.append ("")
        
        for thread in self.threads :
            # find the information
            active_jobs = [job for job in thread[JOBS] if (job[DEBUT] <= value) and (value <= job[FIN] or job[FIN] == -1)]
            active_mutexes = [mutex for mutex in thread[STATESWITCH] if (mutex[NEED] <= value) and (value <= mutex[FREE] or mutex[FREE] == -1)]
            active_tags = [tag for tag in thread[TAGS] if (tag[DT] <= value) and (value <= tag[FT] or tag[FT] == -1)]
            self.info_widget.append (f"<u><b>Thread \"{thread[NOMT]}\"</b></u>")
            
            # display the  information
            if active_jobs :
                for job in active_jobs :
                    if job[NOMJ] == FQNAME:
                        if job[FIN] != -1 :
                            val = f"{job[DEBUT]:,}"
                            val = val.replace(",", " ")
                            val1 = f"{job[FIN]:,}"
                            val1 = val1.replace(",", " ")
                            self.info_widget.append (f"▪️ <b>{job[NOMJ]}</b> (from {val} to {val1}) by {job[ERROR]})")
                        else :
                            val = f"{job[DEBUT]:,}"
                            val = val.replace(",", " ")
                            self.info_widget.append (f"▪️ <b>{job[NOMJ]}</b> (from {val} qutting not finished)")
                    else :
                        if job[FIN] != -1 :
                            val = f"{job[DEBUT]:,}"
                            val = val.replace(",", " ")
                            val1 = f"{job[FIN]:,}"
                            val1 = val1.replace(",", " ")
                            self.info_widget.append (f"▪️ Job <b>{job[NOMJ]} and ID {job[IDJ]}</b> {('child of ' + str (job[IDPARENT])) if job[IDPARENT] > -1 else ''} (from {val} to {val1}) {'✅' if job[SUCCESS] else '❌ : ' + job[ERROR]}")
                        else :
                            val = f"{job[DEBUT]:,}"
                            val = val.replace(",", " ")
                            self.info_widget.append (f"▪️ Job <b>{job[NOMJ]} and ID {job[IDJ]}</b> {('child of ' + str (job[IDPARENT])) if job[IDPARENT] > -1 else ''} (from {val} job not finished)")
            else :
                self.info_widget.append ("▪️ No active job")
            
            if active_mutexes :
                for mutex in active_mutexes :
                    if mutex[HAVE] == -1 :
                        val = f"{mutex[NEED]:,}"
                        val = val.replace(",", " ")
                        self.info_widget.append (f"▪️ Mutex <b>{mutex[NAMEM]}</b> was asked at {val} never taken")
                    elif mutex[FREE] == -1 :
                        val = f"{mutex[NEED]:,}"
                        val = val.replace(",", " ")
                        val1 = f"{mutex[HAVE]:,}"
                        val1 = val1.replace(",", " ")
                        self.info_widget.append (f"▪️ Mutex <b>{mutex[NAMEM]}</b> was asked at {val} and taken at {val1} never freed")
                    else :
                        val = f"{mutex[NEED]:,}"
                        val = val.replace(",", " ")
                        val1 = f"{mutex[HAVE]:,}"
                        val1 = val1.replace(",", " ")
                        val2 = f"{mutex[FREE]:,}"
                        val2 = val2.replace(",", " ")
                        self.info_widget.append (f"▪️ Mutex <b>{mutex[NAMEM]}</b> was asked at {val} and taken at {val1} and freed at {val2}")
            else :
                self.info_widget.append ("▪️ No used mutex")
            
            if active_tags :
                for tag in active_tags :
                    if tag[FT] != -1 :
                        val = f"{tag[DT]:,}"
                        val = val.replace(",", " ")
                        val1 = f"{tag[FT]:,}"
                        val1 = val1.replace(",", " ")
                        self.info_widget.append (f"▪️ Tag <b>{tag[NAMET]}</b> (from {val} to {val1})")
                    else :
                        val = f"{tag[DT]:,}"
                        val = val.replace(",", " ")
                        self.info_widget.append (f"▪️ Tag <b>{tag[NAMET]}</b> (from {val} tag not finished)")
            else :
                self.info_widget.append ("▪️ No active tag")
                
            self.info_widget.append ("")
    
    # -- Methods to create and manege the tracking rectangle
    def create_tracking_rectangle (self) -> None :
        '''
        Create the tracking rectangle with draggable edges
        '''
        # create the tracking rectangle with glass effect
        self.tracking_rect = qtW.QGraphicsRectItem (
            self.left_margin,
            self.top_margin + self.thread_h / 4,
            min (self.width_track_rectangle * 5 * self.pixels_per_unit, self.final_time + self.pixels_per_unit * 2),
            self.nb_threads * self.thread_h
        )
        
        # Gradient for glass effect
        gradient = qtG.QLinearGradient (0, 0, 0, self.nb_threads * self.thread_h)
        color = qtG.QColor (IRON_BLUE_COLOR)
        color.setAlpha (220)  # set the transparency
        gradient.setColorAt (0, color)
        color = qtG.QColor (LAVENDER_COLOR)
        color.setAlpha (180)  # set the transparency
        gradient.setColorAt (1, color)
        
        self.tracking_rect.setBrush (qtG.QBrush (gradient))
        
        # border with dashed line
        color = qtG.QColor (DODGER_BLUE_COLOR)
        color.setAlpha (200)  # set the transparency
        pen = qtG.QPen (color)
        pen.setWidth (2)
        pen.setStyle (qtC.Qt.PenStyle.DashLine)
        pen.setDashPattern ([3, 3])
        
        self.tracking_rect.setPen (pen)
        self.tracking_rect.setZValue (TRACKING_RECTANGLE_ZVALUE)
        self.tracking_rect.setVisible (False)
        self.scene.addItem (self.tracking_rect)

        # shadow effect
        shadow = qtW.QGraphicsDropShadowEffect ()
        shadow.setBlurRadius (15)
        color = qtG.QColor (BLACK_COLOR)
        color.setAlpha (100) # set the transparency
        shadow.setColor (color)
        shadow.setOffset (3, 3)
        self.tracking_rect.setGraphicsEffect (shadow)

        # active the mouse event
        self.tracking_rect.setAcceptHoverEvents (True)
    
    def update_tracking_rectangle (self, value : int) -> None :
        '''
        Update the tracking rectangle to stay centered on the triangle slider and within timeline bounds
        
        Args :
            value (int) : the value of the time line
        '''
        center_x = self.left_margin + (value + 1) * self.pixels_per_unit
        rect_width = self.width_track_rectangle * 10 * self.pixels_per_unit
        
        t_start = self.left_margin
        t_stop = self.left_margin + (self.final_time + 2) * self.pixels_per_unit

        # Calculate the new rectangle boundaries
        new_left = max (t_start, center_x - rect_width / 2)
        new_right = min (t_stop, center_x + rect_width / 2)

        # Adjust the rectangle width to fit within the timeline
        adjusted_width = new_right - new_left

        # Restore the original width if fully within bounds
        if new_left > t_start and new_right < new_right :
            adjusted_width = rect_width
            new_left = center_x - rect_width / 2
            new_right = center_x + rect_width / 2

        self.tracking_rect.setRect (new_left, self.tracking_rect.rect ().y (),
                                    adjusted_width, self.tracking_rect.rect ().height ())
        
        # Ensure the rectangle and markers remain visible if the button is checked
        if self.track_button is not None and self.track_button.isChecked () :
            self.apply_zoom_to_tracking_rect ()
            self.tracking_rect.setVisible (True)
    
    def apply_zoom_to_tracking_rect (self) -> None :
        '''
        Apply zoom effect to the tracking rectangle and its children
        '''
        if self.tracking_rect is None or not self.tracking_rect.isVisible () :
            return

        # Remove existing mask if it exists
        if self.track_mask is not None :
            self.scene.removeItem (self.track_mask)
        
        # create a new mask
        center = self.left_margin + (self.slider.value () + 1) * self.pixels_per_unit
        rect = self.tracking_rect.rect ()
        margin = 2
        inner_rect = qtC.QRectF (rect.x () + margin, rect.y (), rect.width () - 2 * margin, rect.height ())
        
        self.track_mask = qtW.QGraphicsRectItem (inner_rect)
        self.track_mask.setBrush (qtG.QBrush (qtC.Qt.GlobalColor.transparent))
        self.track_mask.setPen (qtG.QPen (qtC.Qt.PenStyle.NoPen))
        self.track_mask.setFlag (qtW.QGraphicsItem.GraphicsItemFlag.ItemClipsChildrenToShape)
        self.track_mask.setZValue (TRACKING_MASK_ZVALUE)
        self.scene.addItem (self.track_mask)

        def in_the_track_rect (item : HighlightableRect) -> None :
            if item is not None and rect.intersects (item.rect ()) :
                self.original_visibility[item] = item.isVisible ()
                
                # create a zoom copy of the item
                item_copy = create_zoomed_copy (item)
                item_copy.setParentItem (self.track_mask)
                
                # Position the copy with the marge
                item_copy.setPos (item.pos ())
                transform = qtG.QTransform ()
                transform.translate (center, rect.center ().y ())
                transform.scale (self.zoom_factor, 1)
                transform.translate (-center, -rect.center ().y ())
                item_copy.setTransform (transform)
                
        for item in self.list_job_item :
            in_the_track_rect (item)
        
        for item in self.list_quitting_item :
            in_the_track_rect (item)
        
        for item1, item2 in self.list_mutex_item :
            in_the_track_rect (item1)
            in_the_track_rect (item2)
        
        for item in self.list_tag_item :
            in_the_track_rect (item)
            
        for item in self.list_deadlock_item :
            in_the_track_rect (item)
    
    def update_width_track_rectangle (self, value : int) -> None :
        '''
        Update the width track rectangle and refresh the tracking rectangle
        
        Args:
            value (int): The new width factor value
        '''
        self.width_track_rectangle = value
        
        if self.tracking_rect is not None and self.tracking_rect.isVisible () :
            self.update_tracking_rectangle (self.slider.value ())

    def update_zoom_factor (self, value : int) -> None :
        '''
        Update the zoom factor and refresh the tracking rectangle
        
        Args:
            value (int): The new zoom factor value
        '''
        self.zoom_factor = value
        self.zoom_slider_label.setText (f"Zoom : {value}x")
        
        if self.tracking_rect is not None and self.tracking_rect.isVisible () :
            self.apply_zoom_to_tracking_rect ()
    
    # -- Methods to show and anshow the path of a job and its parents --
    def show_path (self, job_id : int) -> None :
        '''
        Show the path of a job and its parent
        
        Args :
            job_id (int) : the id of the job we want to show the path
        '''
        self.current_job_id = job_id
        
        # Set all jobs to low opacity
        for job in self.list_job_item :
            job.setOpacity (0.3)
            job.default_pen = job.save_default_pen
            job.setPen (job.default_pen) # Update the pen
            
        for arrow in self.list_arrow :
            self.scene.removeItem (arrow)
        self.list_arrow.clear ()

        # Set opacity for the selected job and its relatives
        current_job = self.list_job_item[job_id]
        current_job.default_pen = current_job.clic_pen
        current_job.setPen (current_job.default_pen) # Update the pen
        current_job.setOpacity (1.0)
        ancetre = self.list_job_parent[job_id][1]
        children = self.list_job_parent[job_id][2]
        
        son = current_job
        nb_recu = 0
        while ancetre > -1 and nb_recu < self.nb_parent:
            parent = self.list_job_item[ancetre]
            parent.setOpacity (1.0)
            
            if self.arrow_button.isChecked () :
                start_point = parent.mapToScene (parent.rect ().x () + 1 * parent.rect ().width (),
                                                parent.rect ().y () + 0.5 * parent.rect ().height ())
                end_point = son.mapToScene (son.rect ().x (),
                                            son.rect ().y () + 0.5 * son.rect ().height ())
                arrow = CurvedArrow (start_point, end_point)
                self.scene.addItem (arrow)
                self.list_arrow.append (arrow)
            
            son = parent
            ancetre = self.list_job_parent[ancetre][1]
            nb_recu += 1
        
        def show_children (parent : int, children : list, nb_recu : int = 0) -> None :
            if children == [] or nb_recu >= self.nb_children :
                return
            for id in children :
                current_job = self.list_job_item[id]
                current_job.setOpacity (1.0)
                
                if self.arrow_button.isChecked () :
                    parent_item = self.list_job_item[parent]
                    start_point = parent_item.mapToScene (parent_item.rect ().x () + 1 * parent_item.rect ().width (),
                                                        parent_item.rect ().y () + 0.5 * parent_item.rect ().height ())
                    end_point = current_job.mapToScene (current_job.rect ().x (),
                                                        current_job.rect ().y () + 0.5 * current_job.rect ().height ())
                    arrow = CurvedArrow (start_point, end_point)
                    self.scene.addItem (arrow)
                    self.list_arrow.append (arrow)
                
                show_children (id, self.list_job_parent[id][2], nb_recu + 1)
        
        show_children (job_id, children)
            
        # Set the opacity of the quitting, mutexes and the tags to 0.3
        for quitting in self.list_quitting_item :
            quitting.setOpacity (0.3)
        for mutex in self.list_mutex_item :
            mutex[0].setOpacity (0.3)
            mutex[1].setOpacity (0.3) if mutex[1] != None else None
        for tag in self.list_tag_item :
            tag.setOpacity (0.3)
    
    def reset_visibility (self) -> None :
        '''
        Reset the visibility of all the view and reset the border of all rectangles
        '''
        # set the opacity of all widget at 1.0
        for item in self.scene.items () :
            if isinstance (item, HighlightableRect) :
                item.setOpacity (1.0)
                item.default_pen = item.save_default_pen
                item.setPen (item.default_pen) # Update the pen
            # Remove the arrows
            if isinstance (item, CurvedArrow) :
                self.scene.removeItem (item)
        self.list_arrow.clear ()
    ################################################################
    
    ################################################################
    # -- Methods for the buttons --
    
    # -- Method which is used when the mutex_button is activate
    def toggle_mutex_view (self, checked : bool) -> None :
        '''
        Toggles the visibility and color scheme of mutex-related items in the GUI.
        When enabled each mutex is assigned a unique color pair.
        When disabled all mutex items are displayed in a default mutex color.

        Args :
            - checked (bool) : A flag indicating whether the mutex view is enabled or disabled.
        '''
        self.mutex_button.setText ("All mutexes" if checked else "Mutex view")
        
        def generate_color_pair (indice) :
            color1 = qtG.QColor (LIST_MUTEX_COLOR[indice % NB_MUTEX_COLORS_VALUE])
            color1.setAlpha (180)
            color2 = qtG.QColor (LIST_MUTEX_COLOR[indice % NB_MUTEX_COLORS_VALUE])
            color2.setAlpha (220)
            return color1, color2

        i = 0
        for m_rect, m_line in self.list_mutex_item :
            name = m_rect.data[NAMEM]  

            if name not in self.mutex_colors :
                self.mutex_colors[name] = generate_color_pair (i)
                i += 1

            if checked :
                color1, color2 = self.mutex_colors[name] 
            else :
                if m_line == None or m_rect.data[HAVE] == -1 :
                    color1, color2 = qtG.QColor (MUTEX_RUNNING_COLOR), qtG.QColor (MUTEX_RUNNING_COLOR)
                else :
                    color1, color2 = qtG.QColor (MUTEX_COLOR), qtG.QColor (MUTEX_COLOR)
                color1.setAlpha (100) # set the transparency
                color2.setAlpha (180) # set the transparency
            m_rect.setBrush (qtG.QBrush (color1))
            m_line.setBrush (qtG.QBrush (color2)) if m_line != None else None
        
        if self.tracking_rect is not None and self.tracking_rect.isVisible () :
            self.update_tracking_rectangle (self.slider.value ())
    
    # -- Method which is used when the arrow_button is activate
    def toggle_arrows (self, checked : bool) -> None :
        '''
        toggle the visibility of the arrows
        
        Args:
            checked (bool) : A flag indicating whether the arrows is enabled or disabled.
        '''
        if checked and self.current_job_id is not None :
            self.show_path (self.current_job_id)
    
    # -- Method which used when the track_button is used
    def toggle_tracking_rectangle (self, checked : bool) -> None :
        '''
        toggle the visibility of the tracking rectangle and their buttons
        
        Args:
            checked (bool) : A flag indicating whether the tracking rectangle is toggle or not
        '''
        self.width_track_rectangle_button.setVisible (checked)
        self.zoom_factor_button.setVisible (checked)
        self.tracking_rect.setVisible (checked)
        
        if checked :
            self.apply_zoom_to_tracking_rect ()
        else : 
            if self.track_mask is not None :
                self.track_mask.setVisible (False)
    
    # -- Methods to manage the number of children and parents shown in the path --
    # Parent methods
    def increase_parent_show (self, nb_parent_max : int) -> None :
        '''
        Increase the number of parent who are show in show_path méthode
        
        Args :
            - nb_parent_max (int) : the maximum number of parent shown
        '''
        if self.nb_parent < nb_parent_max :
            self.nb_parent += 1
            self.text_value_parent.setText (str (self.nb_parent))
            if self.current_job_id is not None :
                self.show_path (self.current_job_id)
    
    def decrease_parent_show (self, nb_parent_min : int) -> None :
        '''
        Decrease the number of parent who are show in show_path méthode
        
        Args :
            - nb_parent_min (int) : the minimum number of parent shown
        '''
        if self.nb_parent > nb_parent_min :
            self.nb_parent -= 1
            self.text_value_parent.setText (str (self.nb_parent))
            if self.current_job_id is not None :
                self.show_path (self.current_job_id)
        
    def update_parent_show (self) -> None :
        '''
        Take the text write by the user to use it, this value must be a int positive between 0 and 9999
        '''
        text = self.text_value_parent.text ()
        if text :
            self.nb_parent = int (text)
            if self.current_job_id is not None :
                self.show_path (self.current_job_id)
                
    # Children methods
    def increase_child_show (self, nb_children_max : int) -> None :
        '''
        Increase the number of child who are show in show_path méthode
        
        Args :
            - nb_children_max (int) : the maximum number of children shown
        '''
        if self.nb_children < nb_children_max :
            self.nb_children += 1
            self.text_value_children.setText (str (self.nb_children))
            if self.current_job_id is not None :
                self.show_path (self.current_job_id)

    def decrease_child_show (self, nb_children_min : int) -> None :
        '''
        Decrease the number of child who are show in show_path méthode

        Args :
            - nb_children_min (int) : the minimum number of children shown
        '''
        if self.nb_children > nb_children_min :
            self.nb_children -= 1
            self.text_value_children.setText (str (self.nb_children))
            if self.current_job_id is not None :
                self.show_path (self.current_job_id)
        
    def update_child_show (self) -> None :
        '''
        Take the text write by the user to use it, this value must be a int positive between 0 and 9999
        '''
        text = self.text_value_children.text ()
        if text :
            self.nb_children = int (text)
            if self.current_job_id is not None :
                self.show_path (self.current_job_id)
    ################################################################
    
    ################################################################
    # -- Methods to display on the scene --
    
    # -- Display the time line --
    def display_time_line (self) -> None :
        '''
        Display the time line on the scene
        '''
        for i in range (self.final_time + 3) :
            pos_x = self.left_margin + i * self.pixels_per_unit
            color = qtG.QColor (GRID_LINE_COLOR)
            
            if i != 0 and i != self.final_time + 2 and (i - 1) % (NB_SHOW_LABEL_TIME_LINE_VALUE * 2) == 0:
                line = self.scene.addLine (pos_x, self.top_margin, pos_x, self.top_margin + (self.nb_threads + 0.5) * self.thread_h, color)
                line.setZValue (GRID_LINE_ZVALUE)
            
            if i <= self.final_time and i % NB_SHOW_LABEL_TIME_LINE_VALUE == 0:
                value = f"{i:,}"
                value = value.replace(",", " ")
                label_t = qtW.QGraphicsSimpleTextItem (value)
                pos_label = pos_x + self.pixels_per_unit - label_t.boundingRect ().width () / 2
                label_t.setPos (pos_label, self.top_margin - 40)
                self.scene.addItem (label_t)
    
    # -- Display the threads row
    def display_threads_row (self) -> None :
        '''
        Display the threads row on the scene
        This method creates a horizontal line for each thread, alternating colors for better visibility.
        '''
        for i in range (self.nb_threads + 2) :
            pos_y = self.top_margin + i * self.thread_h - ((self.thread_h * 3 / 4) if i != 0 else 0)
            bg_color = qtG.QColor (ROW_THREAD_0_COLOR) if i % 2 == 0 else qtG.QColor (ROW_THREAD_1_COLOR)
            
            bg_rect = qtW.QGraphicsRectItem (self.left_margin, pos_y, (self.final_time + 2) * self.pixels_per_unit, self.thread_h if (i != (self.nb_threads + 2) - 1 and i != 0) else self.thread_h / 4)
            bg_rect.setBrush (bg_color)
            bg_rect.setPen (qtG.QPen (qtC.Qt.PenStyle.NoPen))
            bg_rect.setZValue (THREADS_ROW_ZVALUE)
            self.scene.addItem (bg_rect)
            
            if i != 0 :
                line = self.scene.addLine (self.left_margin, pos_y, self.left_margin + (self.final_time + 2) * self.pixels_per_unit, pos_y, qtG.QColor (GRID_LINE_COLOR))
                line.setZValue (GRID_LINE_ZVALUE)
    
    # -- Display the threads --
    def display_threads (self) -> None :
        '''
        Display the threads on the scene
        '''   
        for i, thread in enumerate (self.threads) : 
            pos_y = i * self.thread_h + self.top_margin + self.thread_h / 4
            
            thread_label = qtW.QGraphicsSimpleTextItem (thread[NOMT])
            thread_label.setPos (self.left_margin - 90, pos_y + (self.thread_h - thread_label.boundingRect ().height ()) / 2)
            thread_label.setZValue (THREADS_LABEL_ZVALUE)
            thread_label.setFlag (qtW.QGraphicsItem.GraphicsItemFlag.ItemIsMovable, False)
            thread_label.setFlag (qtW.QGraphicsItem.GraphicsItemFlag.ItemIsSelectable, False)
            self.scene.addItem (thread_label)
            
            # Display the jobs
            self.display_jobs (thread, pos_y)
            # Display the mutexes
            self.display_mutexes (thread, pos_y)
            # Display the tags
            self.display_tags (thread, pos_y)
    
    # -- Display the jobs --
    def display_jobs (self, thread : dict, thread_pos : float) -> None :
        '''
        Display the jobs on the scene
        
        Args :
            - thread (dict) : the thread data containing the jobs
            - thread_pos (float) : the position of the thread in the scene
        '''
        for j, job in enumerate (thread[JOBS]) :
            # Add information, if is not a quitting job, in the tree of jobs
            if job[NOMJ] != FQNAME :
                self.list_job_parent[job[IDJ]][0] = thread[IDT]
                self.list_job_parent[job[IDJ]][1] = job[IDPARENT] 
                (self.list_job_parent[job[IDPARENT]][2].append (job[IDJ])) if job[IDPARENT] > -1 else None
            
            if job[FIN] == -1 : # job running
                fin = self.final_time
            else : # Normal job
                fin = job[FIN]
            
            j_x = self.pixels_per_unit + self.left_margin + job[DEBUT] * self.pixels_per_unit
            width = (fin - job[DEBUT]) * self.pixels_per_unit + (0 if (fin - job[DEBUT] != 0) else 1) # to avoid a width of 0
            
            if job[NOMJ] == FQNAME : # Quitting
                j_y = thread_pos
                height = self.thread_h
                color = qtG.QColor(JOB_QUITTING_COLOR)
                pen = qtG.QPen (qtG.QColor (JOB_QUITTING_COLOR), 3)
                hpen = qtG.QPen (qtG.QColor (JOB_HIGHLIGHT_COLOR), 4)
                
                j_rect = HighlightableRectQuitting (j_x, j_y, width, height,
                                            data = job, 
                                            hover_label = self.hover_label,
                                            default_pen = pen, 
                                            highlight_pen = hpen)
                j_rect.setZValue (QUITTING_ZVALUE)
            
            elif job[FIN] == -1 : # job running
                j_y = thread_pos + 5
                height = self.thread_h - 10
                color = qtG.QColor (JOB_RUNNING_COLOR)
                pen = qtG.QPen (qtG.QColor (JOB_RUNNING_COLOR), 1)
                hpen = qtG.QPen (qtG.QColor (JOB_HIGHLIGHT_COLOR), 2)
                
                j_rect = HighlightableRectJob (j_x, j_y, width, height, 
                                            parent_window = self.window, 
                                            data = job, 
                                            hover_label = self.hover_label, 
                                            callback = lambda id_j = job[IDJ] : self.show_path (id_j),
                                            default_pen = pen, 
                                            highlight_pen = hpen)
                j_rect.setZValue (JOB_ZVALUE)
                
            elif job[DEBUT] == job[FIN] : # job instant
                j_y = thread_pos + 4.5
                height = self.thread_h - 9
                color = qtG.QColor (JOB_INSTANT_OK_COLOR) if job[SUCCESS] else qtG.QColor (JOB_INSTANT_KO_COLOR) # OK job (vert) else KO job (rouge)
                pen = qtG.QPen (qtG.QColor (JOB_INSTANT_OK_COLOR), 4) if job[SUCCESS] else qtG.QPen (qtG.QColor (JOB_INSTANT_KO_COLOR), 4) 
                hpen = qtG.QPen (qtG.QColor (JOB_HIGHLIGHT_COLOR), 5)
                
                j_rect = HighlightableRectJob (j_x, j_y, width, height, 
                                            parent_window = self.window, 
                                            data = job, 
                                            hover_label = self.hover_label, 
                                            callback = lambda id_j = job[IDJ] : self.show_path (id_j),
                                            default_pen = pen, 
                                            highlight_pen = hpen)
                j_rect.setZValue (JOB_INSTANT_ZVALUE)
                
            else : # job non instant
                j_y = thread_pos + 5
                height = self.thread_h - 10
                color = qtG.QColor (JOB_OK_COLOR) if job[SUCCESS] else qtG.QColor (JOB_KO_COLOR)  # OK job (vert) else KO job (rouge)
                pen = qtG.QPen (qtG.QColor (JOB_BORDER_COLOR), 1)
                hpen = qtG.QPen (qtG.QColor (JOB_HIGHLIGHT_COLOR), 2)
                
                j_rect = HighlightableRectJob (j_x, j_y, width, height, 
                                            parent_window = self.window, 
                                            data = job, 
                                            hover_label = self.hover_label, 
                                            callback = lambda id_j = job[IDJ] : self.show_path (id_j),
                                            default_pen = pen, 
                                            highlight_pen = hpen)
                j_rect.setZValue (JOB_ZVALUE)
            
            j_rect.setBrush (qtG.QBrush (color))
            j_rect.setPen (pen)
            j_rect.setFlag (qtW.QGraphicsItem.GraphicsItemFlag.ItemIsMovable, False)
            
            self.scene.addItem (j_rect)
            if job[NOMJ] == FQNAME :
                self.list_quitting_item.append (j_rect)
            else :
                self.list_job_item[job[IDJ]] = j_rect
    
    # -- Display the mutexes --
    def display_mutexes (self, thread : dict, thread_pos : float) -> None :
        '''
        Display the mutexes on the scene
        
        Args :
            - thread (dict) : the thread data containing the mutexes
            - thread_pos (float) : the position of the thread in the scene
        '''
        for k, mutex in enumerate (thread[STATESWITCH]) :
            # Need -> Have area
            if mutex[HAVE] == -1 : # Mutex running
                fin = self.final_time
                color = qtG.QColor (MUTEX_RUNNING_COLOR)
            else :
                fin = mutex[HAVE] # Normal mutex
                color = qtG.QColor (MUTEX_COLOR)

            m_nh_x = self.pixels_per_unit + self.left_margin + mutex[NEED] * self.pixels_per_unit
            width = (fin - mutex[NEED]) * self.pixels_per_unit + (0 if (fin - mutex[NEED] != 0) else 1) # to avoid 0 width

            m_nh_y = thread_pos + self.thread_h / 5
            height = self.thread_h - (self.thread_h / 5) * 2
            
            hpen = qtG.QPen (color, 1)
            color.setAlpha (100)
            pen = qtG.QPen (color, 0)
            
            m_nh_rect = HighlightableRectMutex (m_nh_x, m_nh_y, width, height, 
                                                data = mutex, 
                                                hover_label = self.hover_label, 
                                                default_pen = pen, 
                                                highlight_pen = hpen)
            m_nh_rect.setBrush (qtG.QBrush (color))
            m_nh_rect.setPen (pen)
            
            m_nh_rect.setZValue (MUTEXES_ZVALUE)
            
            m_nh_rect.setFlag (qtW.QGraphicsItem.GraphicsItemFlag.ItemIsMovable, False)
            self.scene.addItem (m_nh_rect)
            
            # Have -> Free area
            if mutex[HAVE] == -1 : # Mutex running
                m_hf_rect = None
            else :
                if mutex[FREE] == -1 : # Mutex running
                    fin = self.final_time
                    color = qtG.QColor (MUTEX_RUNNING_COLOR)
                else : # normal mutex
                    fin = mutex[FREE]
                    color = qtG.QColor (MUTEX_COLOR)
                    
                m_hf_x = self.pixels_per_unit + self.left_margin + mutex[HAVE] * self.pixels_per_unit
                width = (fin - mutex[HAVE]) * self.pixels_per_unit + (0 if (fin - mutex[HAVE] != 0) else 1) # to avoid 0 width
                
                m_hf_y = thread_pos + self.thread_h / 6
                height = self.thread_h - (self.thread_h / 6) * 2
                
                hpen = qtG.QPen (color, 1)
                color.setAlpha (180)
                pen = qtG.QPen (color, 0)
                
                m_hf_rect = HighlightableRectMutex (m_hf_x, m_hf_y, width, height, 
                                            data = mutex, 
                                            hover_label = self.hover_label, 
                                            default_pen = pen, 
                                            highlight_pen = hpen)
                m_hf_rect.setBrush (qtG.QBrush (color))
                m_hf_rect.setPen (pen)
                
                m_hf_rect.setZValue (MUTEXES_ZVALUE)
                
                m_hf_rect.setFlag (qtW.QGraphicsItem.GraphicsItemFlag.ItemIsMovable, False)
                self.scene.addItem (m_hf_rect)
            
            self.list_mutex_item.append ((m_nh_rect, m_hf_rect))
    
    # -- Display the tags --
    def display_tags (self, thread : dict, thread_pos : float) -> None :
        '''
        Display the tags on the scene
        
        Args :
            - thread (dict) : the thread data containing the tags
            - thread_pos (float) : the position of the thread in the scene
        '''
        for t, tag in enumerate (thread[TAGS]) :
            if tag[FT] == -1 : # Tag running
                fin = self.final_time
                color = qtG.QColor (TAG_RUNNING_COLOR)
                pen = qtG.QPen (qtG.QColor (TAG_RUNNING_COLOR), 1)
            else : # Normal tag
                fin = tag[FT]
                color = qtG.QColor (TAG_COLOR)
                pen = qtG.QPen (qtG.QColor (TAG_COLOR), 1)
            
            t_x = self.pixels_per_unit + self.left_margin + tag[DT] * self.pixels_per_unit
            t_width = (fin - tag[DT]) * self.pixels_per_unit + (0 if (fin - tag[DT] != 0) else 1) # to avoid 0 width
            
            t_y = thread_pos + self.thread_h / 2 - 2.5
            t_height = 5
            
            hpen = qtG.QPen (qtG.QColor (TAG_HIGHLIGHT_COLOR), 2)
            
            t_rect = HighlightableRectTag (t_x, t_y, t_width, t_height, 
                                            data = tag, 
                                            hover_label = self.hover_label,
                                            default_pen = pen, 
                                            highlight_pen = hpen)
            t_rect.setBrush (qtG.QBrush (color))
            t_rect.setPen (pen)
            
            t_rect.setZValue (TAGS_ZVALUE)
            t_rect.setFlag (qtW.QGraphicsItem.GraphicsItemFlag.ItemIsMovable, False)
            self.scene.addItem (t_rect)
            self.list_tag_item.append (t_rect)
            
    # -- Display the deadlock situation --
    def display_deadlock (self) -> None :
        '''
        Display the deadlock situations on the scene
        '''
        for th in self.deadlock[DEADLOCK_CYCLE] :
            for i, thread in enumerate (self.threads) :
                if thread[IDT] == th[DEADLOCK_TID] :
                    thread_name = thread[NOMT]
                    thread_pos = i
                d_x = self.left_margin + (self.deadlock[DEADLOCK_TIME] + 1) * self.pixels_per_unit - 1
                d_y = i * self.thread_h + self.top_margin + self.thread_h / 4
                d_width = 2
                d_height = self.thread_h
                color = qtG.QColor (DEADLOCK_COLOR)
                pen = qtG.QPen (qtG.QColor (DEADLOCK_COLOR), 1)
                hpen = qtG.QPen (qtG.QColor (DEADLOCK_HIGHLIGHT_COLOR), 2)
                
                d_rect = HighlightableRectDeadlock (d_x, d_y, d_width, d_height,
                                                    data = self.data,
                                                    hover_label = self.hover_label,
                                                    default_pen = pen, 
                                                    highlight_pen = hpen)

                d_rect.setBrush (qtG.QBrush (color))
                d_rect.setPen (pen)
                
                d_rect.setZValue (DEADLOCK_ZVALUE)
                d_rect.setFlag (qtW.QGraphicsItem.GraphicsItemFlag.ItemIsMovable, False)
                self.scene.addItem (d_rect)
                self.list_deadlock_item.append (d_rect)
    ################################################################

def main () -> None :
    # Force new instance
    app = qt.QtWidgets.QApplication.instance ()
    if app is None :
        app = qt.QtWidgets.QApplication (sys.argv)
    
    # Clear existing window
    for widget in qt.QtWidgets.QApplication.allWidgets () :
        widget.deleteLater ()
    
    # Create new window
    window = MainWindow ()
    window.show ()
    app.exec ()

if __name__ == "__main__" :
    main ()

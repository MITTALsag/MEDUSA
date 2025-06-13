################################################################
# Constant for the GUI
################################################################

################################################################
# -- name of the keys in the json file --
FINAL_TIME = 'final_time'
NB_THREADS = 'nb_threads'
NB_JOBS_OK = 'nb_jobs_OK'
NB_JOBS_KO = 'nb_jobs_KO'
NB_JOBS_RUNNING = 'nb_running_jobs'
NB_QUITTING = 'nb_quitting'
THREADS = 'threads'
DEADLOCK = 'DeadLock'

NOMT = 'nomT'
IDT = 'idT'
JOBS = 'jobs'
STATESWITCH = 'StateSwitch'
TAGS = 'Tags'

NOMJ = 'nomJ'
IDJ = 'id'
IDPARENT = 'id parent'
DEBUT = 'debut'
FIN = 'fin'
SUCCESS = 'success'
ERROR = 'error'

NAMEM = 'NameM'
NEED = 'NEED'
HAVE = 'HAVE'
FREE = 'FREE'

NAMET = 'NameT'
DT = 'START'
FT = 'STOP'

DEADLOCK_TIME = 'Time'
DEADLOCK_CYCLE = 'Cycle'
DEADLOCK_TID = 'ThreadId'
DEADLOCK_NEED = 'Need'
DEADLOCK_HAVE = 'Have'

# specific name in the json file
FQNAME = 'Quitting'
################################################################

################################################################
# -- Colors for the GUI -- hexadecimal format Red Green Blue -- 
WHITE_COLOR = 0xFFFFFF # white
BLACK_COLOR = 0x000000 # black
IRON_BLUE_COLOR = 0xDCE6FA # iron blue
LAVENDER_COLOR = 0xE6DCFA # lavender
DODGER_BLUE_COLOR = 0x1E90FF # dodger blue
TRANSPARENT_COLOR = "transparent" # if you want change to an another color, you need to change the code

DEFAULT_HIGHTLIGHT_COLOR = 0x0000FF # navy blue
DEFAULT_PEN_COLOR = 0x000000 # black
DEFAULT_BRUSH_COLOR = 0x000000 # black

MENU_BACKGROUND_COLOR = 0xFFDD99 # light yellow
DATA_BACKGROUND_COLOR = 0x99D6FF # light blue

BUTTON_BACKGROUND_COLOR = 0xADD8E6 # lightblue
BUTTON_HOVER_COLOR = 0x0070B9 # blue
BUTTON_CHECKED_COLOR = 0x0070B9 # blue
BUTTON_HOVER_TEXT_COLOR = 0xFFFFFF # white
BUTTON_CHECKED_TEXT_COLOR = 0xFFFFFF # white
BUTTON_BORDER_COLOR = 0x0070B9 # blue

TOOLTIP_BACKGROUND_COLOR = 0xFFFFFF # white
TOOLTIP_BORDER_COLOR = 0x000000 # black

GRID_LINE_COLOR = 0xD3D3D3 # light gray
ROW_THREAD_0_COLOR = 0xE0E9EB # light gray
ROW_THREAD_1_COLOR = 0xFFFFFF # white

JOB_HIGHLIGHT_COLOR = 0x000080 # navy blue
JOB_QUITTING_COLOR = 0xFFA500 # orange
JOB_RUNNING_COLOR = 0xa000a0 # purple
JOB_INSTANT_OK_COLOR = 0x1CCD19 # green
JOB_INSTANT_KO_COLOR = 0xCD1919 # red
JOB_OK_COLOR = 0x8BDF77 # light green
JOB_KO_COLOR = 0xDF7777 # light red
JOB_BORDER_COLOR = 0x000000 # black

MUTEX_COLOR = 0x7F8071 # gray
MUTEX_RUNNING_COLOR = 0x939485 # dark gray

TAG_HIGHLIGHT_COLOR = 0x000080 # navy blue
TAG_COLOR = 0x0000FF  # blue
TAG_RUNNING_COLOR = 0x600060 # dark blue

DEADLOCK_COLOR = 0xA52A2A # red brown
DEADLOCK_HIGHLIGHT_COLOR = 0x76381A # brown

LIST_MUTEX_COLOR = [
    0xFF1744, # vivid red
    0x2979FF, # vivid blue
    0x00E676, # vivid green
    0xFFD600, # vivid yellow
    0xFF9100, # vivid orange
    0xD500F9, # vivid purple
    0x00B8D4, # vivid cyan
    0xFF4081, # vivid pink
    0x76FF03, # vivid lime
    0xC51162, # vivid magenta
    0x00BFAE, # vivid teal
    0xFFAB00, # vivid amber
    0x304FFE, # vivid indigo
    0x64DD17, # vivid light green
    0xAEEA00, # vivid chartreuse
    0xFF6D00, # vivid deep orange
    0x6200EA, # vivid deep purple
    0x0091EA, # vivid light blue
    0x00C853, # vivid emerald
    0xFF5252, # vivid coral
] # Modifie the NB_MUTEX_COLORS_VALUE si tu modifies le nombre de couleurs
################################################################

################################################################
# -- Icons names .png --
FILE_ICON = "icon_gui"
GENIUS_EFFECT_ICON = "loupe.png"
FOUR_DIRECTIONNEL_ARROW_ICON = "fleche-de-direction.png"
TWO_DIRECTIONNEL_ARROW_ICON = "bidirectionnel.png"
RESET_ZOOM_ICON = "zoom.png"
################################################################

################################################################
# -- Default values --
TOP_MARGIN_VALUE = 50
LEFT_MARGIN_VALUE = 100 # apply that name of threads are short
PIXELS_PER_UNIT_VALUE = 5
THREADS_HEIGHT_MIN_VALUE = 50
NB_SHOW_LABEL_TIME_LINE_VALUE = 50

INIT_WIDTH_TRACK_RECTAGLE_VALUE = 5
MIN_WIDTH_TRACK_RECTAGLE_VALUE = 1
MAX_WIDTH_TRACK_RECTAGLE_VALUE = 20
INIT_ZOOM_VALUE = 5
MIN_ZOOM_VALUE = 2
MAX_ZOOM_VALUE = 20
INIT_CHILD_PATH_VALUE = 2
INIT_PARENT_PATH_VALUE = 2
MIN_PATH_VALUE = 0
MAX_PATH_VALUE = 9999

NB_MUTEX_COLORS_VALUE = 20
################################################################

################################################################
# -- Set Z values --
THREADS_ROW_ZVALUE = 0
GRID_LINE_ZVALUE = 1

THREADS_LABEL_ZVALUE = 2
JOB_ZVALUE = 2
QUITTING_ZVALUE = 3
JOB_INSTANT_ZVALUE = 4
MUTEXES_ZVALUE = 5
TAGS_ZVALUE = 6
DEADLOCK_ZVALUE = 7

ARROW_ZVALUE = 8

TRIANGLE_SLIDER_ZVALUE = 10
TRACKING_RECTANGLE_ZVALUE = 11
SLIDER_LINE_ZVALUE = 10
TRACKING_MASK_ZVALUE = 20
################################################################

################################################################
# -- Emojis --
JOB_OK_EMOJI = '🟩'
JOB_KO_EMOJI = '🟥'
JOB_RUNNING_EMOJI = '🟪'
MUTEX_EMOJI = '⬛'
TAG_EMOJI = '🟦'
QUITTING_EMOJI = '🟧'
DEADLOCK_EMOJI = '🟫'
################################################################

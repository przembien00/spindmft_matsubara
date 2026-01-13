black = 0x000000
white = 0xFFFFFF
grey = 0x4C4C4C
darkgrey = 0x3B3B3B
darkblue = 0x0C283D
blue = 0x0000FF
orange = 0xFF8800
lightorange = 0xFFD900
red = 0xFF0000
cfblue = 0x6495ED
cbblue = 0x0047AB
lightblue = 0x7EB6E0

fau_blue        = 0x002F6C
fau_darkblue    = 0x041E42
phil_yellow     = 0xFFB81C
phil_orange     = 0xE87722
rw_red          = 0xC8102E
rw_darkred      = 0x971B2F
med_blue        = 0x00A3E0
med_darkblue    = 0x0061A0
nat_green       = 0x43B02A
nat_darkgren    = 0x228848
tf_metallic     = 0x779FB5
tf_darkmetallic = 0x41748D

tuverylightgreen = 0xE3F6BB
tulightgreen = 0xBDEA61
tugreen     = 0x83B818
tudarkgreen = 0x669013
tuverydarkgreen = 0x50710F
tuorange    = 0xFFAF00
tublue      = 0x4169E1 # royalblue
tdo_orangedark  = 0xCA7406

def hex_to_rgb( hex_value ):
    b = (hex_value & 0xFF) / 255.0
    g = ((hex_value >> 8) & 0xFF) / 255.0
    r = ((hex_value >> 16) & 0xFF) / 255.0
    return r, g, b

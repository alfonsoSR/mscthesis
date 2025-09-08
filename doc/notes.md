# Questions Dominic

- Definition of shape, rotation, and shape-deformation settings
- My target for the pre-fit residuals is something close to the order of magnitude of the expected noise I extract from the IFMS/FDETs. This is something I'd like to discuss.
- Access to marcopolo: comment conversation with Giuseppe and implications for first steps of phase A
- Discussion of patterns on the closed-loop pre-fit residuals produced by Luigi's code,


# Updates

- Update condition from A to A' considering that the degree to which open and cl match might not be good enough for B.

- Plot patterns in residual over longer period of time to try to identify the origin of the patterns. Use estrack stations if working with IFMS files


# misc

- Chat with michael about time ephemerides
- Could have a look at impact of lorentz transformation of spacecraft coordinates


# SPICE kernels

- CK is about attitude of spacecraft structures or instruments

EARTHSTNS_ITRF93_YYMMDD.BSP    SPICE Kernel (SPK) that contains ephemeris
                                   data for NASA DSN stations relative to the
                                   terrestrial reference frame label 'ITR93'.
                                   This file was released on YY-MM-DD.

ESTRACK_Vvv.BSP                SPICE Kernel (SPK) that contains the
                                   position for each of the ESA ESTRACK ground
                                   stations relative to the center of the
                                   Earth. In the interest of flexibility, in
                                   this file the reference frame is labeled
                                   with the alias 'EARTH_FIXED'. Any
                                   application using this file must map the
                                   alias 'EARTH_FIXED' to either 'ITRF93' or
                                   'IAU_EARTH'.

MEX_ROB_YYMMDD_yymmdd_vvv.BSP  SPICE Kernel (SPK) that contains high
                                   precision ephemeris for Mars Express (MEX),
                                   created by the Royal Observatory of
                                   Belgium. The kernels cover MEX position
                                   from YY-MM-DD to yy-mm-dd; vvv is the
                                   version number.

MEX_STRUCT_Vvv.BSP             SPICE Kernel (SPK) that contains the
                                   location of various Mars-Express
                                   instruments and structures.

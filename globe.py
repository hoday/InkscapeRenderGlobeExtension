#!/usr/bin/env python
'''
Globe rendering extension for Inkscape
Copyright (C) 2009 Gerrit Karius

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


About the Globe rendering extension:



'''
from __future__ import division
import inkex
import simplestyle, sys
from math import *

def draw_SVG_line(x1, y1, x2, y2, width, name, parent):
    style = { 'stroke': '#000000', 'stroke-width':str(width), 'fill': 'none' }
    line_attribs = {'style':simplestyle.formatStyle(style),
                    inkex.addNS('label','inkscape'):name,
                    'd':'M '+str(x1)+','+str(y1)+' L '+str(x2)+','+str(y2)}
    inkex.etree.SubElement(parent, inkex.addNS('path','svg'), line_attribs )
    
def draw_SVG_rect(x,y,w,h, width, fill, name, parent):
    style = { 'stroke': '#000000', 'stroke-width':str(width), 'fill':fill}
    rect_attribs = {'style':simplestyle.formatStyle(style),
                    inkex.addNS('label','inkscape'):name,
                    'x':str(x), 'y':str(y), 'width':str(w), 'height':str(h)}
    inkex.etree.SubElement(parent, inkex.addNS('rect','svg'), rect_attribs )

def draw_SODIPODI_elipse(cx,cy,rx,ry, width, fill, name, parent):
    style = { 'stroke': '#000000', 'stroke-width':str(width), 'fill':fill}
    circle_attribs = {'style':simplestyle.formatStyle(style),
                      inkex.addNS('label','inkscape'):name,
                      inkex.addNS('cx','sodipodi'):str(cx),
                      inkex.addNS('cy','sodipodi'):str(cy),
                      inkex.addNS('rx','sodipodi'):str(rx),
                      inkex.addNS('ry','sodipodi'):str(ry),
                      inkex.addNS('type','sodipodi'):'arc'}
    inkex.etree.SubElement(parent, inkex.addNS('path','svg'), circle_attribs)    

def draw_SODIPODI_elipse_rotated(cx,cy,rx,ry, width, fill, name, parent, rotationAngle):
    a = cos(rotationAngle)
    b = sin(rotationAngle)
    c = -sin(rotationAngle)
    d = cos(rotationAngle)
    e = -(a*cx + c*cy) + cx
    f = -(b*cx + d*cy) + cy
    style = { 'stroke': '#000000', 'stroke-width':str(width), 'fill':fill}
    if rx == 0:
        x1 = cx 
        x2 = cx
        y1 = cy - ry
        y2 = cy + ry
        circle_attribs = {'style':simplestyle.formatStyle(style),
                          inkex.addNS('label','inkscape'):name,
                          'd':'M '+str(x1)+','+str(y1)+' L '+str(x2)+','+str(y2),
                          'transform':'matrix('+str(a)+','+str(b)+','+str(c)+','+str(d)+','+str(e)+','+str(f)+')'}
    elif ry == 0: 
        x1 = cx - rx 
        x2 = cx + rx
        y1 = cy
        y2 = cy
        circle_attribs = {'style':simplestyle.formatStyle(style),
                          inkex.addNS('label','inkscape'):name,
                          'd':'M '+str(x1)+','+str(y1)+' L '+str(x2)+','+str(y2),
                          'transform':'matrix('+str(a)+','+str(b)+','+str(c)+','+str(d)+','+str(e)+','+str(f)+')'}                           
    else:
        circle_attribs = {'style':simplestyle.formatStyle(style),
                          inkex.addNS('label','inkscape'):name,
                          inkex.addNS('cx','sodipodi'):str(cx),
                          inkex.addNS('cy','sodipodi'):str(cy),
                          inkex.addNS('rx','sodipodi'):str(rx),
                          inkex.addNS('ry','sodipodi'):str(ry),
                          inkex.addNS('type','sodipodi'):'arc',
                          'transform':'matrix('+str(a)+','+str(b)+','+str(c)+','+str(d)+','+str(e)+','+str(f)+')'}
    inkex.etree.SubElement(parent, inkex.addNS('path','svg'), circle_attribs)

def draw_SODIPODI_elipse_segment_rotated(cx,cy,rx,ry, width, fill, name, parent, rotationAngle, segmentAngleStart, segmentAngleEnd):
    a = cos(rotationAngle)
    b = sin(rotationAngle)
    c = -sin(rotationAngle)
    d = cos(rotationAngle)
    e = -(a*cx + c*cy) + cx
    f = -(b*cx + d*cy) + cy
    style = { 'stroke': '#000000', 'stroke-width':str(width), 'fill':fill}
    if rx == 0:
        x1 = cx 
        x2 = cx
        y1 = cy - ry
        y2 = cy + ry
        circle_attribs = {'style':simplestyle.formatStyle(style),
                          inkex.addNS('label','inkscape'):name,
                          'd':'M '+str(x1)+','+str(y1)+' L '+str(x2)+','+str(y2),
                          'transform':'matrix('+str(a)+','+str(b)+','+str(c)+','+str(d)+','+str(e)+','+str(f)+')'}
    elif ry == 0: 
        x1 = cx - rx 
        x2 = cx + rx
        y1 = cy
        y2 = cy
        circle_attribs = {'style':simplestyle.formatStyle(style),
                          inkex.addNS('label','inkscape'):name,
                          'd':'M '+str(x1)+','+str(y1)+' L '+str(x2)+','+str(y2),
                          'transform':'matrix('+str(a)+','+str(b)+','+str(c)+','+str(d)+','+str(e)+','+str(f)+')'}     
    else:  
        circle_attribs = {'style':simplestyle.formatStyle(style),
                          inkex.addNS('label','inkscape'):name,
                          inkex.addNS('cx','sodipodi'):str(cx),
                          inkex.addNS('cy','sodipodi'):str(cy),
                          inkex.addNS('rx','sodipodi'):str(rx),
                          inkex.addNS('ry','sodipodi'):str(ry),
                          inkex.addNS('start','sodipodi'):str(segmentAngleStart),
                          inkex.addNS('end','sodipodi'):str(segmentAngleEnd),
                          inkex.addNS('open','sodipodi'):'true',
                          inkex.addNS('type','sodipodi'):'arc',
                          'transform':'matrix('+str(a)+','+str(b)+','+str(c)+','+str(d)+','+str(e)+','+str(f)+')'}
    inkex.etree.SubElement(parent, inkex.addNS('path','svg'), circle_attribs)


class Globe(inkex.Effect):
    def __init__(self):
        inkex.Effect.__init__(self)
        self.OptionParser.add_option("--longitudeLineCount",
                        action="store", type="int", 
                        dest="longitudeLineCount", default=8,
                        help="Number of longitude lines")
        self.OptionParser.add_option("--latitudeLineCount",
                        action="store", type="int", 
                        dest="latitudeLineCount", default=5,
                        help="Number of latitude lines")
        self.OptionParser.add_option("--rotationXDegrees",
                        action="store", type="float", 
                        dest="rotationXDegrees", default=0,
                        help="Rotation around X axis (degrees)")
        self.OptionParser.add_option("--rotationYDegrees",
                        action="store", type="float", 
                        dest="rotationYDegrees", default=0,
                        help="Rotation around Y axis (degrees)")
        self.OptionParser.add_option("--isSeeThrough",
                        action="store", type="inkbool", 
                        dest="isSeeThrough", default=False,
                        help="Is the globe see-through")
    def effect(self):

        name = 'globe'

        # globe fill and stroke style
        fill = 'none'
        width = 1

        #input parameters - globe center and radius
        cyb = 500.0
        cxb = 500.0 
        rb  = 100.0 

        longitudeRotationAngleDegrees = float(self.options.rotationYDegrees)
        tiltForwardAngleDegrees       = float(self.options.rotationXDegrees)

        # inputs range fixing
        # tiltForwardAngle is adjusted to vary from 0 <= angle < pi
        if tiltForwardAngleDegrees >= 180.0:
            tiltForwardAngleDegrees -= 180.0
        elif tiltForwardAngleDegrees < 180.0:
            tiltForwardAngleDegrees += 180.0 



        if self.options.longitudeLineCount > 0:
            angleSpacingLongitudeLinesDegrees = 180.0 / float(self.options.longitudeLineCount);
            # longitudeAngle is wrapped to vary from 0 <= angle < angleSpacingLongitudeLines.
            while longitudeRotationAngleDegrees < 0:
                longitudeRotationAngleDegrees += angleSpacingLongitudeLinesDegrees
            while longitudeRotationAngleDegrees >= angleSpacingLongitudeLinesDegrees:
                longitudeRotationAngleDegrees -= angleSpacingLongitudeLinesDegrees                   

        # units conversion from degrees to radians
        tiltForwardAngle = tiltForwardAngleDegrees * pi / 180.0;
        initialAngleLongitudeLines = longitudeRotationAngleDegrees * pi / 180.0

        # derived parameters
        rxb = rb
        ryb = rb

        #
        # start drawing
        #

        # create the group to put the globe in
        group_attribs = {inkex.addNS('label','inkscape'):name}
        parent = inkex.etree.SubElement(self.current_layer, 'g', group_attribs)


        # draw the outside border
        draw_SODIPODI_elipse_rotated(cxb,cyb,rxb,ryb, width, fill, 'border', parent, 0)

        # draw the longitude lines
        # elipse #0 corresponds to ring on the front (visible only as a straight vertical line)
        # elipse #n-1 corresponds to the ring that is almost 180 degrees away
        # elipse #n/2 corresponds to ring around the side (overlaps with globe boundary) (only if n is even)
        if self.options.longitudeLineCount > 0:
            angleSpacingLongitudeLines = pi / float(self.options.longitudeLineCount);
            yOfPole = ryb * cos(tiltForwardAngle)
            for i in range(0, self.options.longitudeLineCount):
                lineName = 'longitude' + str(i)
                # longitudeAngle is always from 0 to pi.
                # rotation angle is always from 0 to pi.
                # rx is never negative.
                longitudeAngle = ((float(i)) * angleSpacingLongitudeLines) + initialAngleLongitudeLines
                if tiltForwardAngleDegrees == 0 or tiltForwardAngleDegrees == 180.0:
                    if longitudeAngle < pi/2:
                        rotationAngle = 0.0
                    else:
                        rotationAngle = pi
                    rx = rxb * sin(longitudeAngle)

                    arcStart = pi/2
                    arcEnd = -pi/2

                else:
                    rotationAngle = acos(cos(longitudeAngle) / sqrt(1 - pow(sin(longitudeAngle)*cos(tiltForwardAngle), 2)))
                    rx = rxb  * sin(longitudeAngle) * cos(tiltForwardAngle)
                    if rx < 0:
                        rx = -rx
                        arcStart = -pi/2
                        arcEnd = pi/2
                    else:
                        arcStart = pi/2
                        arcEnd = -pi/2
                ry = ryb
                cx = cxb
                cy = cyb
                if self.options.isSeeThrough:
                    draw_SODIPODI_elipse_rotated(cx,cy,rx,ry, width, fill, lineName, parent, rotationAngle)
                else:
                    draw_SODIPODI_elipse_segment_rotated(cx,cy,rx,ry, width, fill, lineName, parent, rotationAngle, arcStart, arcEnd)


        # draw the latitude lines
        # elipse #0 corresponds to ring closest to north pole.
        # elipse #n-1 corresponds to ring closest to south pole.
        # equator is ring #(n-1)/2 (only if n is odd). 
        if self.options.latitudeLineCount > 0:
            angleSpacingLatitudeLines  = pi / (1.0 + float(self.options.latitudeLineCount));
            yOfPole = ryb * cos(tiltForwardAngle)
            for i in range(0, self.options.latitudeLineCount):
                lineName = 'latitude' + str(i)
                # angleOfCurrentLatitudeLine is always from 0 to pi.
                # tiltForwardAngle is always from 0 to pi.
                # ry is never negative.                
                angleOfCurrentLatitudeLine = float(i + 1) * angleSpacingLatitudeLines
                rx = rxb * sin(angleOfCurrentLatitudeLine)
                ry = rx * sin(tiltForwardAngle)
                cx = cxb
                cy = cyb - yOfPole*cos(angleOfCurrentLatitudeLine)
                if self.options.isSeeThrough:
                    draw_SODIPODI_elipse_rotated(cx,cy,rx,ry, width, fill, lineName, parent, 0)
                else:
                    if tiltForwardAngle > pi/2:
                        # tilt away from viewaer
                        if rxb * cos(angleOfCurrentLatitudeLine) / cos(tiltForwardAngle) > rxb:
                            # elipse is not visible
                            pass
                        else:
                            if rxb * cos(angleOfCurrentLatitudeLine) / cos(tiltForwardAngle) < -rxb:
                                # elipse is all visible
                                segmentAngle = pi
                            else:
                                # elipse is only partially visible
                                segmentAngle = acos(max(-1,min(1, -tan(tiltForwardAngle) / tan(angleOfCurrentLatitudeLine))))
                            draw_SODIPODI_elipse_segment_rotated(cx,cy,rx,ry, width, fill, lineName, parent, 0, pi/2+segmentAngle, pi/2-segmentAngle)
                    else:
                        # tilt towards viewer
                        if rxb * cos(angleOfCurrentLatitudeLine) / cos(tiltForwardAngle) < -rxb:
                            # elipse is not visible
                            pass
                        else:
                            if rxb * cos(angleOfCurrentLatitudeLine) / cos(tiltForwardAngle) > rxb:
                                # elipse is all visible
                                segmentAngle = pi
                            else:
                                # elipse is only partially visible
                                segmentAngle = acos(max(-1,min(1, tan(tiltForwardAngle) / tan(angleOfCurrentLatitudeLine))))
                            draw_SODIPODI_elipse_segment_rotated(cx,cy,rx,ry, width, fill, lineName, parent, 0, -pi/2+segmentAngle, -pi/2-segmentAngle)


    def unused(self):
        
        #find the pixel dimensions of the overall grid
        ymax = self.options.dy * self.options.y_divs
        xmax = self.options.dx * self.options.x_divs
        
        # Embed grid in group
        #Put in in the centre of the current view
        t = 'translate(' + str( self.view_center[0]- xmax/2.0) + ',' + \
                           str( self.view_center[1]- ymax/2.0) + ')'
        g_attribs = {inkex.addNS('label','inkscape'):'Grid_Polar:X' + \
                     str( self.options.x_divs )+':Y'+str( self.options.y_divs ),
                     'transform':t }
        grid = inkex.etree.SubElement(self.current_layer, 'g', g_attribs)
        
        #Group for major x gridlines
        g_attribs = {inkex.addNS('label','inkscape'):'MajorXGridlines'}
        majglx = inkex.etree.SubElement(grid, 'g', g_attribs)

        #Group for major y gridlines
        g_attribs = {inkex.addNS('label','inkscape'):'MajorYGridlines'}
        majgly = inkex.etree.SubElement(grid, 'g', g_attribs)
        
        #Group for minor x gridlines
        if self.options.x_subdivs > 1:#if there are any minor x gridlines
            g_attribs = {inkex.addNS('label','inkscape'):'MinorXGridlines'}
            minglx = inkex.etree.SubElement(grid, 'g', g_attribs)
        
        #Group for subminor x gridlines
        if self.options.x_subsubdivs > 1:#if there are any minor minor x gridlines
            g_attribs = {inkex.addNS('label','inkscape'):'SubMinorXGridlines'}
            mminglx = inkex.etree.SubElement(grid, 'g', g_attribs)
        
        #Group for minor y gridlines
        if self.options.y_subdivs > 1:#if there are any minor y gridlines
            g_attribs = {inkex.addNS('label','inkscape'):'MinorYGridlines'}
            mingly = inkex.etree.SubElement(grid, 'g', g_attribs)
        
        #Group for subminor y gridlines
        if self.options.y_subsubdivs > 1:#if there are any minor minor x gridlines
            g_attribs = {inkex.addNS('label','inkscape'):'SubMinorYGridlines'}
            mmingly = inkex.etree.SubElement(grid, 'g', g_attribs)

            
        draw_SVG_rect(0, 0, xmax, ymax, self.options.border_th,
                      'none', 'Border', grid) #border rectangle
        
        #DO THE X DIVISONS======================================
        sd  = self.options.x_subdivs #sub divs per div
        ssd = self.options.x_subsubdivs #subsubdivs per subdiv
        
        for i in range(0, self.options.x_divs): #Major x divisons
            if i>0: #dont draw first line (we made a proper border)
                draw_SVG_line(self.options.dx*i, 0,
                              self.options.dx*i,ymax,
                              self.options.x_divs_th,
                              'MajorXDiv'+str(i), majglx)
            
            if self.options.x_log: #log x subdivs
                for j in range (1, sd):
                    if j>1: #the first loop is only for subsubdivs
                        draw_SVG_line(self.options.dx*(i+log(j, sd)), 0,
                                      self.options.dx*(i+log(j, sd)), ymax,
                                      self.options.x_subdivs_th,
                                      'MinorXDiv'+str(i)+':'+str(j), minglx)
                                  
                    for k in range (1, ssd): #subsub divs
                        if (j <= self.options.x_half_freq) or (k%2 == 0):#only draw half the subsubdivs past the half-freq point
                            if (ssd%2 > 0) and (j > self.options.y_half_freq): #half frequency won't work with odd numbers of subsubdivs,
                                ssd2 = ssd+1 #make even
                            else:
                                ssd2 = ssd #no change
                            draw_SVG_line(self.options.dx*(i+log(j+k/float(ssd2),sd )), 0,
                                          self.options.dx*(i+log(j+k/float(ssd2),sd )), ymax,
                                          self.options.x_subsubdivs_th,'SubminorXDiv'+str(i)+':'+str(j)+':'+str(k), mminglx)
            
            else: #linear x subdivs
                for j in range (0, sd):
                    if j>0: #not for the first loop (this loop is for the subsubdivs before the first subdiv)
                        draw_SVG_line(self.options.dx*(i+j/float(sd)), 0,
                                      self.options.dx*(i+j/float(sd)), ymax,
                                      self.options.x_subdivs_th,
                                      'MinorXDiv'+str(i)+':'+str(j), minglx)
                    
                    for k in range (1, ssd): #subsub divs
                        draw_SVG_line(self.options.dx*(i+(j*ssd+k)/((float(sd)*ssd))) , 0,
                                      self.options.dx*(i+(j*ssd+k)/((float(sd)*ssd))) , ymax,
                                      self.options.x_subsubdivs_th,
                                      'SubminorXDiv'+str(i)+':'+str(j)+':'+str(k), mminglx)
         
        #DO THE Y DIVISONS========================================
        sd  = self.options.y_subdivs    #sub divs per div
        ssd = self.options.y_subsubdivs #subsubdivs per subdiv
                                      
        for i in range(0, self.options.y_divs): #Major y divisons
            if i>0:#dont draw first line (we will make a border)
                draw_SVG_line(0, self.options.dy*i,
                              xmax, self.options.dy*i,
                              self.options.y_divs_th,
                              'MajorYDiv'+str(i), majgly)
            
            if self.options.y_log: #log y subdivs
                for j in range (1, sd):
                    if j>1: #the first loop is only for subsubdivs
                        draw_SVG_line(0,    self.options.dy*(i+1-log(j,sd)),
                                      xmax, self.options.dy*(i+1-log(j,sd)),
                                      self.options.y_subdivs_th,
                                      'MinorXDiv'+str(i)+':'+str(j), mingly)
                    
                    for k in range (1, ssd): #subsub divs
                        if (j <= self.options.y_half_freq) or (k%2 == 0):#only draw half the subsubdivs past the half-freq point
                            if (ssd%2 > 0) and (j > self.options.y_half_freq): #half frequency won't work with odd numbers of subsubdivs,
                                ssd2 = ssd+1
                            else:
                                ssd2 = ssd #no change
                            draw_SVG_line(0,    self.options.dx*(i+1-log(j+k/float(ssd2),sd )),
                                          xmax, self.options.dx*(i+1-log(j+k/float(ssd2),sd )),
                                          self.options.y_subsubdivs_th,
                                          'SubminorXDiv'+str(i)+':'+str(j)+':'+str(k), mmingly)
            else: #linear y subdivs
                for j in range (0, self.options.y_subdivs):
                    if j>0:#not for the first loop (this loop is for the subsubdivs before the first subdiv)
                        draw_SVG_line(0,    self.options.dy*(i+j/float(sd)),
                                      xmax, self.options.dy*(i+j/float(sd)),
                                      self.options.y_subdivs_th,
                                      'MinorXYiv'+str(i)+':'+str(j), mingly)
                    
                    for k in range (1, ssd): #subsub divs
                        draw_SVG_line(0,    self.options.dy*(i+(j*ssd+k)/((float(sd)*ssd))),
                                      xmax, self.options.dy*(i+(j*ssd+k)/((float(sd)*ssd))),
                                      self.options.y_subsubdivs_th,
                                      'SubminorXDiv'+str(i)+':'+str(j)+':'+str(k), mmingly)



if __name__ == '__main__':
    e = Globe()
    e.affect()



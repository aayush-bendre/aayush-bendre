import arcade

class MyWindow(arcade.Window):
    def __init__(self, width, height, windowtitle):
        super().__init__(width, height, windowtitle)
        self.x = width/2
        self.y = height/2
        self.dx = 1
        self.dy = 1
        self.color = (0,255,255)
        

        arcade.set_background_color( (0,0,0) )

    def on_draw(self):
        arcade.start_render()
        arcade.draw_rectangle_filled(self.x, self.y, 64, 64, self.color)

   
    def on_update(self, dt):
        self.x += 100*self.dx*dt
        self.y += 100*self.dy*dt

        if self.x > self.width - 32:
            self.x = self.width - 32
            self.dx = -self.dx
        
        if self.x < 32:
            self.x = 32
            self.dx = -self.dx

        if self.y > self.height - 32:
            self.y = self.height - 32
            self.dy = -self.dy
        
        if self.y < 32:
            self.y = 32
            self.dy = -self.dy
        

    
window = MyWindow(640,480,"Title")
arcade.run()

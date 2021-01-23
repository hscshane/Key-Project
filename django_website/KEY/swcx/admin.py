from django.contrib import admin
from .models import SWCX 

@admin.register(SWCX)
class KEYAdmin(admin.ModelAdmin):
    pass 

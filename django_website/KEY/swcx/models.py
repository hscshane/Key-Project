from django.db import models

class SWCX(models.Model):
    name = models.CharField(max_length=100)
    description = models.TextField(blank=True)

